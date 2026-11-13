#!/usr/bin/env python3
"""
Comprehensive statistical analysis script for CapScreen data.

This script performs QC and DESeq2 analysis on CapScreen count data.

1. Variant metrics, following the screening notebooks:
   - Sum technical replicates within each (group, bio_rep) before QC and DESeq2
   - RPM and log2(RPM + 1) from the full library, then a per-group low-count filter
     (a group is low when every replicate is < min_count, including zeros)
   - CV of RPM across biological replicates

2. QC plots and tables, written to ``output_dir/qc``:
   - log2(RPM) Pearson scatter plots and a replicate-correlation table
   - CV density histograms and an AAV2 RPM contamination bar plot
   - Blind VST sample-sample Pearson heatmap
   - PCA on blind VST after dropping the Input group (all positive-variance
     variants, and the top variable variants)

3. DESeq2, written to ``output_dir/deseq2``:
   - Contrasts are taken only from Sample_info ``enrichment_role`` /
     ``enrichment_arm`` (target versus each background on an arm)
   - Each contrast is fit pairwise on only the samples in those two groups
     (preferred when groups differ in variability, e.g. low-variance Input)
   - Design is ``~ group``, or ``~ batch_seq + group`` when multiple
     sequencing batches are present and not confounded with group
   - Wald statistics, log2 fold changes, and BH-adjusted p-values

4. Enrichment / enhancement (Lin et al.), written to
   ``output_dir/enrichment_enhancement``:
   - Same Sample_info role/arm columns as DESeq2
   - Er = target RPM / input RPM; Eh = target RPM / background RPM
   - Median- and mean-based scores, hit calling, Er-Eh scatter plots
   - Top enriched variants: grouped log10 RPM bar plots for all groups
"""

import argparse
import json
import logging
import re
import sys
from pathlib import Path
from typing import Tuple, Optional, Dict, Any, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from itertools import combinations
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import LogFormatterMathtext, LogLocator, NullLocator
from scipy.stats import pearsonr
from sklearn.decomposition import PCA

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Set up logging
logger = logging.getLogger(__name__)


def load_config(config_file: Optional[Path] = None) -> Dict[str, Any]:
    """
    Load configuration from JSON file.
    
    Args:
        config_file: Path to config.json file. If None, looks for config.json in capscreen directory.
        
    Returns:
        Dictionary with statistical analysis configuration values
    """
    if config_file is None:
        # Try to find config.json in capscreen directory
        capscreen_dir = Path(__file__).parent.parent
        config_file = capscreen_dir / "config.json"
        if not config_file.exists():
            config_file = Path.cwd() / "config.json"
    
    if not config_file.exists():
        logger.warning(f"Config file not found at {config_file}, using default values")
        # Return default config
        return {
            "analysis": {
                "n_cpus": 8,
                "low_expression_threshold": 10,
                "pca": {
                    "top_n_features": 500
                }
            },
            "differential_expression": {
                "min_replicates": 2,
                "min_count": 10,
                "drop_if_any_group_low": False,
                "log2rpm_pseudocount": 1,
                "adjust_batch": True,
            },
            "enrichment_enhancement": {
                "er_cut": 3,
                "eh_cut": 3,
                "min_rpm": 5,
                "n_top": 10,
                "pseudocount": 1,
            },
            "output": {
                "write_excel": True,
            },
            "plots": {
                "dpi": {
                    "standard": 300
                },
                "figure_sizes": {
                    "pca": [7.5, 6.5]
                },
                "pca": {
                    "point_size": 120,
                    "fontsize": 10
                }
            }
        }
    
    with open(config_file, 'r') as f:
        full_config = json.load(f)
    
    # Extract statistical_analysis section, or use full config if section doesn't exist
    config = full_config.get("statistical_analysis", full_config)
    
    logger.info(f"Loaded configuration from {config_file}")
    return config


def setup_logging(verbose: bool = False, log_file: Optional[Path] = None):
    """Set up logging to stdout and optionally a dedicated stat log file."""
    level = logging.DEBUG if verbose else logging.INFO
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    root = logging.getLogger()
    root.setLevel(level)
    has_stream = any(
        isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
        for h in root.handlers
    )
    if not has_stream:
        stream = logging.StreamHandler(sys.stdout)
        stream.setFormatter(fmt)
        stream.setLevel(level)
        root.addHandler(stream)
    if log_file:
        _attach_stat_log_file(root, Path(log_file), fmt, level)


def _attach_stat_log_file(
    logger_obj: logging.Logger,
    log_file: Path,
    fmt: Optional[logging.Formatter] = None,
    level: int = logging.INFO,
) -> None:
    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    abs_path = str(log_file.resolve())
    for handler in logger_obj.handlers:
        if isinstance(handler, logging.FileHandler) and getattr(handler, "baseFilename", None) == abs_path:
            return
    if fmt is None:
        fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(fmt)
    file_handler.setLevel(level)
    logger_obj.addHandler(file_handler)


def _stat_failure_reason(step: str, exc: BaseException) -> str:
    """One-line reason for a failed STAT step."""
    name = type(exc).__name__
    msg = str(exc).strip() or "(no error message)"
    if isinstance(exc, MemoryError) or "Unable to allocate" in msg:
        return (
            f"{name}: {msg}. Out of memory during '{step}'. "
            "Increase Docker/host RAM (8-16 GB recommended) and re-run capscreen stat."
        )
    return f"{name}: {msg}"


def read_count_table(count_file: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read merged count table and separate raw counts from RPM counts.
    
    Args:
        count_file: Path to merged.counts.csv file
        
    Returns:
        Tuple of (raw_counts_df, rpm_counts_df)
    """
    logger.info(f"Reading count table from {count_file}")
    
    # Read the count table
    counts = pd.read_csv(count_file, sep=',')
    if "ID" not in counts.columns and "ID_WLG" in counts.columns:
        counts = counts.rename(columns={"ID_WLG": "ID"})
    
    # Filter out unassigned peptides — they are not valid library members and
    # would create duplicate index entries that break downstream analysis.
    if 'ID' in counts.columns:
        n_before = len(counts)
        counts = counts[~counts['ID'].str.lower().eq('unassigned')]
        n_removed = n_before - len(counts)
        if n_removed > 0:
            logger.info(f"Excluded {n_removed:,} unassigned variant(s) from statistical analysis")
    
    # Set index to ID
    if 'ID' in counts.columns:
        counts.index = counts['ID']
    else:
        logger.warning("ID column not found, using first column as index")
        counts.index = counts.iloc[:, 0]
    
    # Identify RPM columns
    rpm_columns = [col for col in counts.columns if col.endswith('_RPM')]
    
    # Create raw counts dataframe (drop RPM columns and metadata columns)
    raw_counts = counts.drop(columns=['ID', 'peptide', 'variable_seq'] + rpm_columns, errors='ignore')
    
    # Create RPM counts dataframe (keep only RPM columns)
    rpm_counts = counts[['ID', 'peptide'] + rpm_columns].copy() if rpm_columns else pd.DataFrame()
    
    # Remove RPM suffix from column names for RPM dataframe
    if not rpm_counts.empty:
        rpm_counts.columns = [col.replace('_RPM', '') if col.endswith('_RPM') else col 
                             for col in rpm_counts.columns]
    
    # Transpose raw counts (samples as rows, variants as columns)
    raw_counts = raw_counts.T
    
    # Remove first row if it's AAV2 or other control (common pattern)
    if raw_counts.index[0] == 'AAV2' or 'AAV' in str(raw_counts.index[0]):
        logger.info("Removing AAV2 control row")
        raw_counts = raw_counts.iloc[1:]
    
    logger.info(f"Raw counts shape: {raw_counts.shape}")
    logger.info(f"RPM counts shape: {rpm_counts.shape}")
    
    return raw_counts, rpm_counts


def read_sample_info(sample_info_file: Path) -> pd.DataFrame:
    """
    Read sample information CSV file.
    
    Args:
        sample_info_file: Path to Sample_info.csv file
        
    Returns:
        DataFrame with sample metadata
    """
    logger.info(f"Reading sample info from {sample_info_file}")
    
    meta = pd.read_csv(sample_info_file, sep=",")
    if "sample_id" not in meta.columns:
        meta = meta.rename(columns={meta.columns[0]: "sample_id"})
    meta["sample_id"] = meta["sample_id"].astype(str)
    if meta["sample_id"].duplicated().any():
        n_dup = int(meta["sample_id"].duplicated().sum())
        logger.warning("Dropping %d duplicate sample_id row(s); keeping the first", n_dup)
        meta = meta.drop_duplicates("sample_id", keep="first")

    required_cols = ["sample_id", "group", "tech_rep", "bio_rep", "batch_seq"]
    missing_cols = [col for col in required_cols if col not in meta.columns]
    if missing_cols:
        logger.warning(f"Missing columns in sample info: {missing_cols}")
    if "group" not in meta.columns or "bio_rep" not in meta.columns:
        raise ValueError("sample_info must contain group and bio_rep columns")
    missing_group = meta["group"].isna() | meta["group"].astype(str).str.strip().eq("")
    if missing_group.any():
        logger.warning("Dropping %d sample_info row(s) with an empty group", int(missing_group.sum()))
        meta = meta.loc[~missing_group].copy()
    if meta.empty:
        raise ValueError("sample_info has no rows with a group label")

    optional_cols = ["enrichment_role", "enrichment_arm"]
    keep = [col for col in required_cols + optional_cols if col in meta.columns]
    meta = meta.loc[:, keep].copy()
    if "enrichment_role" in meta.columns:
        meta["enrichment_role"] = meta["enrichment_role"].fillna("").astype(str).str.strip()
    if "enrichment_arm" in meta.columns:
        meta["enrichment_arm"] = meta["enrichment_arm"].fillna("").astype(str).str.strip()
    meta = meta.set_index("sample_id", drop=False)
    meta.index.name = "sample_id"
    
    logger.info(f"Sample info shape: {meta.shape}")
    logger.info(f"Groups found: {sorted(meta['group'].astype(str).unique())}")
    return meta


def identify_groups(meta: pd.DataFrame, config: Optional[Dict[str, Any]] = None) -> Tuple[list, list, list]:
    """
    Identify target, background, and input groups from metadata.
    
    First checks for explicit group role definitions in config file.
    Falls back to string heuristics if config is not provided or incomplete.
    
    Config format (in statistical_analysis.group_roles):
    {
        "input": ["Input_group1", "Input_group2"],
        "target": ["Target_group1", "CD47_WT"],
        "background": ["Bead_group1", "CD47_KO"]
    }
    
    If group_roles is empty or not provided, uses heuristics:
    - Input: groups containing "input" (case-insensitive)
    - Target: groups containing "target", "wt", "cd", "receptor"
    - Background: groups containing "bead", "ko", "background"
    
    Args:
        meta: DataFrame with sample metadata
        config: Optional configuration dictionary with group role definitions.
                Should contain "group_roles" key with "input", "target", "background" lists.
        
    Returns:
        Tuple of (target_groups, background_groups, input_groups)
    """
    all_groups = meta['group'].unique().tolist()
    
    # Check for explicit group role definitions in config
    group_roles = None
    if config is not None:
        group_roles = config.get("group_roles", None)
    
    # Check if group_roles exists and has at least one non-empty list
    # If all lists are empty, fall back to heuristics
    has_explicit_roles = (group_roles is not None and 
                         isinstance(group_roles, dict) and
                         any(isinstance(v, list) and len(v) > 0 
                             for v in group_roles.values()))
    
    if has_explicit_roles:
        # Use explicit configuration
        input_groups = group_roles.get("input", [])
        target_groups = group_roles.get("target", [])
        background_groups = group_roles.get("background", [])
        
        # Validate that configured groups exist in metadata
        input_groups = [g for g in input_groups if g in all_groups]
        target_groups = [g for g in target_groups if g in all_groups]
        background_groups = [g for g in background_groups if g in all_groups]
        
        # Warn about missing groups
        missing_input = [g for g in group_roles.get("input", []) if g not in all_groups]
        missing_target = [g for g in group_roles.get("target", []) if g not in all_groups]
        missing_background = [g for g in group_roles.get("background", []) if g not in all_groups]
        
        if missing_input:
            logger.warning(f"Configured input groups not found in metadata: {missing_input}")
        if missing_target:
            logger.warning(f"Configured target groups not found in metadata: {missing_target}")
        if missing_background:
            logger.warning(f"Configured background groups not found in metadata: {missing_background}")
        
        # Find unassigned groups (groups in metadata but not in config)
        assigned_groups = set(input_groups + target_groups + background_groups)
        unassigned_groups = [g for g in all_groups if g not in assigned_groups]
        
        if unassigned_groups:
            logger.warning(f"Groups found in metadata but not in config group_roles: {unassigned_groups}")
            logger.info("Unassigned groups will be classified using heuristics")
            # Classify unassigned groups using heuristics
            for g in unassigned_groups:
                role = _classify_group_heuristic(g)
                if role == "input":
                    input_groups.append(g)
                elif role == "target":
                    target_groups.append(g)
                else:
                    background_groups.append(g)
        
        logger.info(f"Using explicit group role configuration from config file")
        logger.info(f"Input groups: {input_groups}")
        logger.info(f"Target groups: {target_groups}")
        logger.info(f"Background groups: {background_groups}")
        
        return target_groups, background_groups, input_groups
    
    # Fall back to heuristics if no config provided
    logger.info("No explicit group role configuration found, using heuristics")
    return _identify_groups_heuristic(meta)


def _classify_group_heuristic(group_name: str) -> str:
    """
    Classify a single group using string heuristics.
    
    Args:
        group_name: Name of the group
        
    Returns:
        "input", "target", or "background"
    """
    group_lower = group_name.lower()
    
    # Check for input groups
    if 'input' in group_lower:
        return "input"
    
    # Check for target groups
    target_keywords = ['target', 'wt', 'WT', 'receptor']
    if any(kw in group_lower for kw in target_keywords):
        return "target"
    
    # Check for background groups
    background_keywords = ['bead', 'KO', 'background']
    if any(kw in group_lower for kw in background_keywords):
        return "background"
    
    # Default to background if unclear
    return "background"


def _identify_groups_heuristic(meta: pd.DataFrame) -> Tuple[list, list, list]:
    """
    Identify target, background, and input groups using string heuristics.
    
    Args:
        meta: DataFrame with sample metadata
        
    Returns:
        Tuple of (target_groups, background_groups, input_groups)
    """
    all_groups = meta['group'].unique().tolist()
    
    # Identify input groups (groups containing "input" or "Input")
    input_groups = [g for g in all_groups if 'input' in g.lower() or 'Input' in g]
    
    # Remaining groups are target/background
    non_input_groups = [g for g in all_groups if g not in input_groups]
    
    # Try to identify target and background groups based on common naming patterns
    # Target patterns: Target, WT, specific receptor names (e.g., CD47)
    # Background patterns: Bead, KO, Bead-Fc
    target_keywords = ['Target', 'target', 'WT', 'wt']
    background_keywords = ['Bead', 'bead', 'KO', 'ko', 'Background', 'background']
    
    target_groups = []
    background_groups = []
    
    for g in non_input_groups:
        is_target = any(kw in g for kw in target_keywords)
        is_background = any(kw in g for kw in background_keywords)
        
        if is_target and not is_background:
            target_groups.append(g)
        elif is_background:
            background_groups.append(g)
        else:
            # If unclear, check if it's a specific receptor name (likely target)
            # or if it contains common background terms
            if any(x in g for x in ['CD47', 'CD', 'receptor']):
                target_groups.append(g)
            else:
                # Default to background if unclear
                background_groups.append(g)
    
    # If we couldn't identify any, put all in target_groups
    if not target_groups and not background_groups:
        target_groups = non_input_groups
    
    logger.info(f"Input groups: {input_groups}")
    logger.info(f"Target groups: {target_groups}")
    logger.info(f"Background groups: {background_groups}")
    
    return target_groups, background_groups, input_groups


def compute_rpm_enrichment_analysis(rpm_counts: pd.DataFrame, meta: pd.DataFrame,
                                    input_groups: list, LAMBDA: float = 0.01,
                                    target_groups: Optional[list] = None,
                                    background_groups: Optional[list] = None,
                                    config: Optional[Dict[str, Any]] = None) -> pd.DataFrame:
    """
    Compute RPM-based enrichment analysis following the notebook pattern.
    
    Steps:
    1. Compute mean RPM across technical replicates per biological replicate and group
    2. Compute log2 enrichment vs Input
    3. Aggregate across experiments
    4. Create summary statistics and target specificity metrics
    
    Args:
        rpm_counts: RPM counts dataframe (variants as rows, samples as columns)
        meta: Metadata dataframe
        input_groups: List of input group names
        LAMBDA: Small constant to avoid divide-by-zero (default: 0.01)
        target_groups: Optional list of target group names (if None, will identify from enrichment columns)
        background_groups: Optional list of background group names (if None, will identify from enrichment columns)
        
    Returns:
        Annotated RPM dataframe with enrichment metrics
    """
    logger.info("Computing RPM-based enrichment analysis...")
    
    # Refactored: Fixed index/column naming to use ID as identifier and preserve peptide column
    # Sourcery-AI suggestion: Keep ID as index (don't rename to 'peptide') and preserve 'peptide' 
    # column (contains actual sequence) to avoid confusion between identifiers and sequences
    # Prepare RPM counts: set ID as index if it's a column, preserve peptide column
    df_rpm = rpm_counts.copy()
    if 'ID' in df_rpm.columns:
        df_rpm = df_rpm.set_index('ID')
    # Keep ID as index name (don't rename to 'peptide')
    # Preserve 'peptide' column (contains actual sequence) - don't drop it
    
    # Remove AAV2 or control variants
    control_variants = df_rpm.index[df_rpm.index.str.contains('AAV', case=False, na=False)]
    if len(control_variants) > 0:
        logger.info(f"Removing control variants from RPM analysis: {list(control_variants)}")
        df_rpm = df_rpm.drop(control_variants)
    
    # Align with metadata
    common_samples = df_rpm.columns.intersection(meta.index)
    df_rpm = df_rpm[common_samples]
    meta_aligned = meta.loc[common_samples]
    
    logger.info(f"Using {len(common_samples)} samples for RPM enrichment analysis")
    
    # Step 1: Compute mean RPM across technical replicates, per biological replicate
    logger.info("Step 1: Computing mean RPM across technical replicates...")
    mu = []
    
    for rep in sorted(meta_aligned["bio_rep"].unique()):
        meta_r = meta_aligned[meta_aligned["bio_rep"] == rep]
        
        for group in meta_aligned["group"].unique():
            samples = meta_r[meta_r["group"] == group].index
            
            if len(samples) == 0:
                continue
            
            # Get samples that exist in df_rpm
            samples_in_rpm = [s for s in samples if s in df_rpm.columns]
            if len(samples_in_rpm) == 0:
                continue
            
            mu_r = df_rpm[samples_in_rpm].mean(axis=1)  # μ_{i,s}
            # Refactored: Changed from "peptide" to "ID" to use identifier consistently
            mu.append(
                pd.DataFrame({
                    "ID": mu_r.index,
                    "mu": mu_r.values,
                    "group": group,
                    "bio_rep": rep
                })
            )
    
    mu_df = pd.concat(mu, ignore_index=True)
    logger.info(f"Computed means for {len(mu_df)} ID-group-bio_rep combinations")
    
    # Step 2: Compute log2 enrichment vs Input
    logger.info("Step 2: Computing log2 enrichment vs Input...")
    enrich = []
    
    # Get non-input groups
    non_input_groups = [g for g in meta_aligned["group"].unique() if g not in input_groups]
    
    for rep in sorted(mu_df["bio_rep"].unique()):
        mu_r = mu_df[mu_df["bio_rep"] == rep]
        
        # Get input means for this biological replicate
        mu_input_df = mu_r[mu_r["group"].isin(input_groups)]
        if len(mu_input_df) == 0:
            logger.warning(f"No input samples found for bio_rep {rep}, skipping")
            continue
        
        # Refactored: Changed groupby and set_index from "peptide" to "ID" for identifier consistency
        # Average across input groups if multiple exist
        mu_input = mu_input_df.groupby("ID")["mu"].mean()
        
        for group in non_input_groups:
            mu_g = (
                mu_r[mu_r["group"] == group]
                .set_index("ID")["mu"]
            )
            
            common = mu_g.index.intersection(mu_input.index)
            
            if len(common) > 0:
                e = np.log2((mu_g.loc[common] + LAMBDA) / (mu_input.loc[common] + LAMBDA))
                
                # Refactored: Changed from "peptide" to "ID" for identifier consistency
                enrich.append(
                    pd.DataFrame({
                        "ID": common,
                        "log2_enrichment": e.values,
                        "group": group,
                        "bio_rep": rep
                    })
                )
    
    if len(enrich) == 0:
        logger.warning("No enrichment data computed. Check input groups.")
        return df_rpm
    
    enrich_df = pd.concat(enrich, ignore_index=True)
    logger.info(f"Computed enrichments for {len(enrich_df)} ID-group-bio_rep combinations")
    
    # Step 3: Aggregate across experiments
    logger.info("Step 3: Aggregating across experiments...")
    # Refactored: Changed groupby from "peptide" to "ID" for identifier consistency
    enrich_summary = (
        enrich_df
        .groupby(["ID", "group"])
        .agg(
            mean_log2_enrichment=("log2_enrichment", "mean"),
            sd_log2_enrichment=("log2_enrichment", "std")
        )
        .reset_index()
    )
    
    # Step 4: Create wide format tables and compute statistics
    logger.info("Step 4: Creating summary statistics...")
    df_rpm_mat = df_rpm.copy()
    
    # Refactored: Changed pivot index from "peptide" to "ID" for identifier consistency
    # Pivot mean and sd
    mean_wide = enrich_summary.pivot(index="ID", columns="group", values="mean_log2_enrichment")
    sd_wide = enrich_summary.pivot(index="ID", columns="group", values="sd_log2_enrichment")
    
    # Rename columns
    # Calculate CV (coefficient of variation) with protection against division by zero
    cv_wide = sd_wide / mean_wide.abs().replace(0, np.nan)
    cv_wide = cv_wide.fillna(0)  # Set CV to 0 where mean is 0
    cv_wide = cv_wide.add_prefix("CV_")
    mean_wide = mean_wide.add_prefix("mean_log2Enrich_")
    sd_wide = sd_wide.add_prefix("sd_log2Enrich_")
    
    # Combine
    peptide_annot = pd.concat([mean_wide, sd_wide, cv_wide], axis=1)
    
    # Compute target specificity for ALL target-background pairs
    # First, get the original group names (before prefix was added)
    # mean_wide columns now have prefix "mean_log2Enrich_", so we need to extract group names
    all_group_cols = mean_wide.columns.tolist()
    all_groups_in_enrich = [col.replace("mean_log2Enrich_", "") for col in all_group_cols]
    
    # Use provided target/background groups if available, otherwise identify from column names
    # Check config first for explicit group role definitions
    config_groups = None
    if config is not None:
        config_groups = config.get("group_roles", None)
    
    if target_groups is None or background_groups is None:
        if config_groups:
            # Use explicit configuration if available
            if target_groups is None:
                target_groups = [g for g in config_groups.get("target", []) if g in all_groups_in_enrich]
            if background_groups is None:
                background_groups = [g for g in config_groups.get("background", []) if g in all_groups_in_enrich]
        
        # Fall back to heuristics for any remaining unclassified groups
        unclassified_groups = [g for g in all_groups_in_enrich 
                              if g not in (target_groups or []) and g not in (background_groups or [])]
        
        if unclassified_groups:
            identified_target_groups = []
            identified_background_groups = []
            
            for group_name in unclassified_groups:
                role = _classify_group_heuristic(group_name)
                if role == "target":
                    identified_target_groups.append(group_name)
                elif role == "background":
                    identified_background_groups.append(group_name)
            
            # Merge with config-based groups
            if target_groups is None:
                target_groups = identified_target_groups
            else:
                target_groups.extend(identified_target_groups)
            
            if background_groups is None:
                background_groups = identified_background_groups
            else:
                background_groups.extend(identified_background_groups)
    
    # Filter to only groups that actually have enrichment data
    target_groups = [g for g in target_groups if g in all_groups_in_enrich]
    background_groups = [g for g in background_groups if g in all_groups_in_enrich]
    
    # Create specificity columns for all target-background pairs
    if target_groups and background_groups:
        logger.info(f"Found {len(target_groups)} target group(s): {target_groups}")
        logger.info(f"Found {len(background_groups)} background group(s): {background_groups}")
        specificity_created = 0
        for target_group in target_groups:
            for background_group in background_groups:
                target_col = f"mean_log2Enrich_{target_group}"
                background_col = f"mean_log2Enrich_{background_group}"
                if target_col in peptide_annot.columns and background_col in peptide_annot.columns:
                    # Create specificity column name (handle special characters like "-")
                    target_name = target_group.replace("-", "").replace("_", "")
                    bg_name = background_group.replace("-", "").replace("_", "")
                    specificity_col_name = f"{target_name}_minus_{bg_name}"
                    peptide_annot[specificity_col_name] = (
                        peptide_annot[target_col] - peptide_annot[background_col]
                    )
                    logger.info(f"Created specificity column: {specificity_col_name} ({target_group} - {background_group})")
                    specificity_created += 1
                else:
                    logger.warning(f"Columns not found for specificity calculation: {target_col} or {background_col}")
                    logger.debug(f"Available columns in peptide_annot: {list(peptide_annot.columns)}")
        if specificity_created == 0:
            logger.warning("No specificity columns were created. Check if target/background enrichment columns exist.")
    else:
        logger.warning(f"Could not identify target/background groups. Target groups: {target_groups}, Background groups: {background_groups}")
        logger.info(f"Available enrichment columns: {list(mean_wide.columns)}")
        logger.info(f"All groups found in enrichment: {[col.replace('mean_log2Enrich_', '') for col in mean_wide.columns]}")
    
    # Refactored: Changed merge key from "peptide" to "ID" for identifier consistency
    # Reset index to merge (peptide_annot index is ID)
    peptide_annot = peptide_annot.reset_index()
    
    # Merge into df_rpm (df_rpm_mat index is ID)
    df_rpm_annot = (
        df_rpm_mat
        .reset_index()  # Reset to make index a column (ID) for merging
        .merge(peptide_annot, on="ID", how="left")
        .set_index("ID")
    )
    
    # Step 5: Add per-replicate enrichments in wide format
    logger.info("Step 5: Adding per-replicate enrichments...")
    # Refactored: Changed pivot_table index from "peptide" to "ID" for identifier consistency
    enrich_wide = (
        enrich_df
        .pivot_table(
            index="ID",
            columns=["group", "bio_rep"],
            values="log2_enrichment"
        )
    )
    
    # Flatten column names
    enrich_wide.columns = [
        f"log2Enrich_{g}_rep{r}" for g, r in enrich_wide.columns
    ]
    
    df_rpm_annot = df_rpm_annot.join(enrich_wide, how="left")
    
    logger.info(f"RPM enrichment analysis completed. Final table shape: {df_rpm_annot.shape}")
    
    return df_rpm_annot


def r_make_names(name: str) -> str:
    """Approximate R ``make.names`` for group labels used in column names."""
    cleaned = re.sub(r"[^A-Za-z0-9_]", ".", str(name))
    if not cleaned or cleaned[0].isdigit():
        cleaned = "X" + cleaned
    return cleaned


def _rep_sort_key(value):
    try:
        return (0, float(value))
    except (TypeError, ValueError):
        return (1, str(value))


def _is_human_ntng1_group(group: str) -> bool:
    return str(group).endswith("-NTNG1") and not str(group).endswith("-Ntng1")


def _is_mouse_ntng1_group(group: str) -> bool:
    return str(group).endswith("-Ntng1")


def sanitize_group_label(group: str) -> str:
    """Case-safe file label. Human NTNG1 and mouse Ntng1 stay distinct."""
    group = str(group)
    if group == "Human-NTNG1":
        return "Human.NTNG1"
    if group == "Mouse-Ntng1":
        return "Mouse.Ntng1"
    if _is_human_ntng1_group(group):
        group = group[: -len("-NTNG1")] + "-Human.NTNG1"
    elif _is_mouse_ntng1_group(group):
        group = group[: -len("-Ntng1")] + "-Mouse.Ntng1"
    group = group.replace("-", ".")
    return re.sub(r"[^A-Za-z0-9._]", "_", group)


def de_result_filename(y_group: str, x_group: str) -> str:
    return f"DE_{sanitize_group_label(y_group)}_vs_{sanitize_group_label(x_group)}.tsv"


def _deseq2_pairs_from_enrichment_arms(arms: pd.DataFrame) -> List[Tuple[str, str]]:
    """
    One DESeq2 contrast per target-versus-background pair in each Er/Eh arm.

    Input is used for enrichment Er only and is not added as a DESeq2 reference here.
    """
    pairs: List[Tuple[str, str]] = []
    for _, arm in arms.iterrows():
        target = str(arm["target"])
        for background in str(arm["backgrounds"]).split(";"):
            background = background.strip()
            if not background or background == target:
                continue
            pair = (target, background)
            if pair not in pairs:
                pairs.append(pair)
    return pairs


def build_enrichment_comparisons(
    sample_info: pd.DataFrame,
    group_col: str = "group",
) -> pd.DataFrame:
    """
    Build y-versus-x DESeq2 contrasts only from Sample_info roles/arms.

    Requires ``enrichment_role`` / ``enrichment_arm``. Each arm yields one
    contrast of target versus each background group. y_group is the numerator
    and x_group is the reference.
    """
    arms = _arms_from_enrichment_metadata(sample_info, group_col=group_col)
    pairs = _deseq2_pairs_from_enrichment_arms(arms)
    if not pairs:
        raise ValueError(
            "No DESeq2 contrasts could be built from Sample_info "
            "(need target and at least one background per enrichment_arm)"
        )
    arm_labels: List[str] = []
    for y_group, x_group in pairs:
        label = ""
        for _, arm in arms.iterrows():
            bgs = [b.strip() for b in str(arm["backgrounds"]).split(";") if b.strip()]
            if str(arm["target"]) == y_group and x_group in bgs:
                label = str(arm["tag"])
                break
        arm_labels.append(label)
    out = pd.DataFrame(pairs, columns=["y_group", "x_group"])
    out["mode"] = "sample_info"
    out["arm"] = arm_labels
    logger.info(
        "Built %d DESeq2 contrast(s) from Sample_info enrichment_role/arm",
        len(out),
    )
    return out


def load_merged_counts(count_file: Path) -> pd.DataFrame:
    """Read merged counts with one row per variant. Unassigned rows are removed."""
    logger.info("Reading count table from %s", count_file)
    counts = pd.read_csv(count_file, sep=",")
    if "ID" not in counts.columns and "ID_WLG" in counts.columns:
        counts = counts.rename(columns={"ID_WLG": "ID"})
    if "ID" not in counts.columns:
        raise ValueError("merged counts must contain an ID column")
    n_before = len(counts)
    counts = counts.loc[~counts["ID"].astype(str).str.lower().eq("unassigned")].copy()
    n_removed = n_before - len(counts)
    if n_removed:
        logger.info("Excluded %d unassigned variant(s)", n_removed)
    if counts["ID"].astype(str).duplicated().any():
        raise ValueError("Duplicate ID values in the merged count table")
    logger.info("Count table shape: %s", counts.shape)
    return counts


def _group_is_low(values: np.ndarray, min_count: float) -> np.ndarray:
    """True when every finite replicate count is strictly below ``min_count``.

    Zeros count as low. A group with no finite values is not marked low.
    """
    finite = np.isfinite(values)
    n_finite = finite.sum(axis=1)
    below = finite & (values < min_count)
    return (n_finite > 0) & (below.sum(axis=1) == n_finite)


def update_variant_metrics(
    df: pd.DataFrame,
    sample_info: pd.DataFrame,
    min_count: float = 10,
    drop_if_any_group_low: bool = False,
    pseudocount: float = 1.0,
    min_rpm_post: Optional[float] = None,
    rpm_post_filter_mode: str = "any",
) -> pd.DataFrame:
    """
    Attach RPM and log2(RPM + pseudocount), then drop low-count variants.

    RPM is computed from the full table before row filtering. A group is low
    when every finite replicate count is strictly below ``min_count``
    (zeros included). With ``drop_if_any_group_low=False``, a variant is
    removed only when every group is low, so all-zero variants are dropped.
    """
    sample_ids = [sid for sid in sample_info["sample_id"].astype(str) if sid in df.columns]
    if not sample_ids:
        raise ValueError("No sample_id from sample_info matched columns in the count table")
    missing = sorted(set(sample_info["sample_id"].astype(str)) - set(sample_ids))
    if missing:
        logger.warning("sample_id values with no count column: %s", missing)

    raw = df[sample_ids].apply(pd.to_numeric, errors="coerce")
    n_na = int(raw.isna().sum().sum())
    if n_na:
        logger.warning("Replacing %d missing count value(s) with 0", n_na)
        raw = raw.fillna(0)
    lib_sizes = raw.sum(axis=0).replace(0, np.nan)
    rpm = raw.div(lib_sizes, axis=1) * 1e6
    log2_rpm = np.log2(rpm + pseudocount)

    info = sample_info.loc[sample_info["sample_id"].astype(str).isin(sample_ids), ["sample_id", "group"]]
    low_flags = []
    for _, sub in info.groupby("group", sort=False):
        cols = [sid for sid in sub["sample_id"].astype(str) if sid in raw.columns]
        low_flags.append(_group_is_low(raw[cols].to_numpy(dtype=float), min_count))
    low_mat = np.column_stack(low_flags)
    if drop_if_any_group_low:
        drop = low_mat.any(axis=1)
    else:
        drop = low_mat.all(axis=1)

    keep = ~drop
    out = df.loc[keep].copy()
    rpm = rpm.loc[keep]
    log2_rpm = log2_rpm.loc[keep]
    for sid in sample_ids:
        out[f"{sid}_RPM"] = rpm[sid].to_numpy()
        out[f"{sid}_log2RPM"] = log2_rpm[sid].to_numpy()

    if min_rpm_post is not None:
        rpm_cols = [f"{sid}_RPM" for sid in sample_ids]
        rpm_values = out[rpm_cols].to_numpy(dtype=float)
        if rpm_post_filter_mode == "all":
            post_keep = np.isfinite(rpm_values).all(axis=1) & (rpm_values >= min_rpm_post).all(axis=1)
        else:
            post_keep = (np.nan_to_num(rpm_values, nan=-np.inf) >= min_rpm_post).any(axis=1)
        out = out.loc[post_keep].copy()

    logger.info(
        "Low-count filter (min_count=%s, drop_if_any_group_low=%s): kept %d/%d variants",
        min_count,
        drop_if_any_group_low,
        len(out),
        len(df),
    )
    return out


def collapse_technical_replicates(
    df: pd.DataFrame,
    sample_info: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Sum raw counts across technical replicates within each (group, bio_rep).

    Returns a count table and sample_info with one sample per biological
    replicate. QC metrics and DESeq2 then use biological replicates as the
    unit of analysis, so technical replicates do not inflate degrees of freedom.
    Precomputed ``*_RPM`` / ``*_log2RPM`` columns are dropped; callers recompute
    them after collapsing.
    """
    sample_ids = [sid for sid in sample_info["sample_id"].astype(str) if sid in df.columns]
    if not sample_ids:
        raise ValueError("No sample_id from sample_info matched columns in the count table")

    info = sample_info.loc[sample_info["sample_id"].astype(str).isin(sample_ids)].copy()
    sizes = info.groupby(["group", "bio_rep"], sort=False).size()
    if int(sizes.max()) <= 1:
        logger.info("No technical replicates to collapse (%d samples)", len(info))
        return df, sample_info

    drop_suffix_cols = [
        c for c in df.columns
        if str(c).endswith("_RPM") or str(c).endswith("_log2RPM")
    ]
    meta_cols = [
        c for c in df.columns
        if c not in sample_ids and c not in drop_suffix_cols
    ]
    out = df.loc[:, meta_cols].copy()
    new_rows = []
    used_ids: set = set()
    n_collapsed = 0
    for (group, bio_rep), sub in info.groupby(["group", "bio_rep"], sort=False):
        cols = [sid for sid in sub["sample_id"].astype(str) if sid in df.columns]
        if not cols:
            continue
        raw = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        if len(cols) == 1:
            new_id = cols[0]
            summed = raw[cols[0]].round().astype(int)
        else:
            base = f"{r_make_names(group)}_bio{bio_rep}"
            new_id = base
            suffix = 2
            while new_id in used_ids or new_id in out.columns:
                new_id = f"{base}_{suffix}"
                suffix += 1
            summed = raw.sum(axis=1).round().astype(int)
            n_collapsed += 1
            logger.info(
                "Collapsed %d technical replicates for group=%s bio_rep=%s into %s",
                len(cols),
                group,
                bio_rep,
                new_id,
            )
        out[new_id] = summed.to_numpy()
        used_ids.add(new_id)
        row = sub.iloc[0].to_dict()
        row["sample_id"] = new_id
        if "tech_rep" in row:
            row["tech_rep"] = 1
        new_rows.append(row)

    new_meta = pd.DataFrame(new_rows)
    new_meta["sample_id"] = new_meta["sample_id"].astype(str)
    new_meta = new_meta.set_index("sample_id", drop=False)
    new_meta.index.name = "sample_id"
    logger.info(
        "Technical-replicate collapse: %d samples -> %d biological replicates (%d multi-tech groups)",
        len(info),
        len(new_meta),
        n_collapsed,
    )
    return out, new_meta


def _samples_per_biorep(sample_info: pd.DataFrame, group: str) -> pd.DataFrame:
    """One sample_info row per biological replicate for a group, sorted by bio_rep.

    After :func:`collapse_technical_replicates`, each bio_rep has a single sample.
    If uncollapsed technical replicates remain, the first row for that bio_rep is used.
    """
    sub = sample_info.loc[sample_info["group"].astype(str) == str(group)].copy()
    reps = sorted(sub["bio_rep"].unique(), key=_rep_sort_key)
    rows = []
    for rep in reps:
        hit = sub.loc[sub["bio_rep"] == rep]
        if len(hit):
            if len(hit) > 1:
                logger.warning(
                    "Group '%s' bio_rep %s still has %d samples after collapse; using the first",
                    group,
                    rep,
                    len(hit),
                )
            rows.append(hit.iloc[0])
    if not rows:
        return sub.iloc[0:0]
    return pd.DataFrame(rows)


def add_cv_rpm_by_group(df: pd.DataFrame, sample_info: pd.DataFrame, prefix: str = "CV_RPM") -> pd.DataFrame:
    """Add one CV column per group. CV is sd(RPM) / mean(RPM) across biological replicates."""
    out = df.copy()
    for group in pd.unique(sample_info["group"].astype(str)):
        chosen = _samples_per_biorep(sample_info, group)
        rpm_cols = [f"{sid}_RPM" for sid in chosen["sample_id"].astype(str) if f"{sid}_RPM" in out.columns]
        if len(rpm_cols) < 2:
            logger.warning("Group '%s': need at least 2 biological replicates for CV; skipping", group)
            continue
        values = out[rpm_cols].to_numpy(dtype=float)
        n_finite = np.isfinite(values).sum(axis=1)
        mean = np.divide(np.nansum(values, axis=1), n_finite, out=np.full(len(values), np.nan), where=n_finite > 0)
        centered = values - mean[:, None]
        sq = np.nansum(np.where(np.isfinite(centered), centered ** 2, 0.0), axis=1)
        std = np.sqrt(np.divide(sq, n_finite - 1, out=np.full(len(values), np.nan), where=n_finite > 1))
        cv = np.divide(std, mean, out=np.full(len(values), np.nan), where=np.isfinite(mean) & (mean != 0))
        cv[(n_finite < 2) | ~np.isfinite(mean) | (mean == 0) | ~np.isfinite(std)] = np.nan
        out[f"{prefix}_{r_make_names(group)}"] = cv
    return out


def calc_pearson_biorep_by_group(df: pd.DataFrame, sample_info: pd.DataFrame) -> pd.DataFrame:
    """Pearson correlation of log2(RPM) between biological replicates, one row per group."""
    rows = []
    pair_names = []
    for group in pd.unique(sample_info["group"].astype(str)):
        chosen = _samples_per_biorep(sample_info, group)
        if len(chosen) < 2:
            continue
        cols = []
        labels = []
        for _, row in chosen.iterrows():
            col = f"{row['sample_id']}_log2RPM"
            if col not in df.columns:
                continue
            cols.append(col)
            labels.append(f"bio{row['bio_rep']}")
        if len(cols) < 2:
            continue
        mat = df[cols].apply(pd.to_numeric, errors="coerce")
        mat.columns = labels
        if len(mat) < 3:
            continue
        row_vals = {"group": group}
        for left, right in combinations(labels, 2):
            name = f"{left}_vs_{right}"
            pair = mat[[left, right]].dropna()
            if len(pair) < 3 or pair[left].nunique() < 2 or pair[right].nunique() < 2:
                row_vals[name] = np.nan
            else:
                row_vals[name] = float(pearsonr(pair[left], pair[right])[0])
            if name not in pair_names:
                pair_names.append(name)
        rows.append(row_vals)
    if not rows:
        logger.warning("No groups with at least 2 biological replicates for log2RPM correlation")
        return pd.DataFrame(columns=["group"])
    out = pd.DataFrame(rows)
    for name in pair_names:
        if name not in out.columns:
            out[name] = np.nan
    return out[["group"] + pair_names]


def plot_log2rpm_pearson_pdf(df: pd.DataFrame, sample_info: pd.DataFrame, out_pdf: Path) -> int:
    """Scatter log2(RPM) for every biological-replicate pair within each group.

    Annotate Pearson r and the number of variants only. Do not show a p-value:
    with n equal to the variant count, the test is almost always significant
    even for weak correlations.
    """
    n_pages = 0
    with PdfPages(out_pdf) as pdf:
        for group in pd.unique(sample_info["group"].astype(str)):
            chosen = _samples_per_biorep(sample_info, group)
            usable = []
            for _, row in chosen.iterrows():
                col = f"{row['sample_id']}_log2RPM"
                if col in df.columns:
                    usable.append((col, row["bio_rep"], row["sample_id"]))
            for (c1, b1, s1), (c2, b2, s2) in combinations(usable, 2):
                pair = df[[c1, c2]].apply(pd.to_numeric, errors="coerce").dropna()
                if len(pair) < 3:
                    continue
                x = pair[c1].to_numpy(dtype=float)
                y = pair[c2].to_numpy(dtype=float)
                r_val = float(pearsonr(x, y)[0])
                fig, ax = plt.subplots(figsize=(7, 6))
                ax.scatter(x, y, s=8, alpha=0.35, c="#2C7FB8", linewidths=0)
                if np.nanstd(x) > 0:
                    slope, intercept = np.polyfit(x, y, 1)
                    xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
                    ax.plot(xs, slope * xs + intercept, color="#D95F0E", lw=0.7)
                ax.set_title(f"{group} | bio{b1} vs bio{b2}")
                ax.set_xlabel(f"{s1} log2(RPM+1)")
                ax.set_ylabel(f"{s2} log2(RPM+1)")
                ax.text(
                    0.02,
                    0.98,
                    f"Pearson r = {r_val:.3f}, n = {len(pair)}",
                    transform=ax.transAxes,
                    va="top",
                    fontsize=9,
                )
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
                n_pages += 1
    if n_pages == 0 and out_pdf.exists():
        out_pdf.unlink()
        logger.warning("No log2RPM replicate pairs to plot")
    else:
        logger.info("Saved %d log2RPM Pearson scatter page(s) to %s", n_pages, out_pdf)
    return n_pages


def plot_cv_histogram_pdf(df: pd.DataFrame, out_pdf: Path, cv_high_threshold: float = 0.5, bins: int = 60) -> int:
    """One density histogram per CV_RPM column."""
    cv_cols = sorted(col for col in df.columns if str(col).startswith("CV_RPM_"))
    if not cv_cols:
        logger.warning("No CV_RPM columns to plot")
        return 0
    n_pages = 0
    with PdfPages(out_pdf) as pdf:
        for col in cv_cols:
            values = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            if len(values) < 2:
                logger.warning("Skipping %s: fewer than 2 finite values", col)
                continue
            n_high = int((values > cv_high_threshold).sum())
            pct = 100.0 * n_high / len(values)
            title = str(col)[len("CV_RPM_"):].replace("_", " ").replace(".", " ")
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.hist(values, bins=bins, density=True, color="#4575b4", edgecolor="white", linewidth=0.25)
            ax.set_title(f"CV of RPM across replicates — {title}")
            ax.set_xlabel("CV (SD / mean RPM)")
            ax.set_ylabel("Density")
            ax.text(
                0.98,
                0.98,
                f"n = {len(values)} | CV > {cv_high_threshold:g}: {n_high} ({pct:.1f}%)",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=9,
            )
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            n_pages += 1
    if n_pages == 0 and out_pdf.exists():
        out_pdf.unlink()
        logger.warning("No CV histograms were written")
    else:
        logger.info("Saved CV histograms to %s", out_pdf)
    return n_pages


def plot_aav2_rpm_contamination(df: pd.DataFrame, sample_info: pd.DataFrame, out_pdf: Path, variant_id: str = "AAV2") -> bool:
    """Bar plot of AAV2 RPM by group, with one point per biological replicate."""
    if "ID" not in df.columns:
        logger.warning("Cannot plot AAV2 RPM: ID column is missing")
        return False
    hits = df.index[df["ID"].astype(str) == variant_id]
    if len(hits) == 0:
        logger.warning("Variant '%s' not found; skipping contamination plot", variant_id)
        return False
    if len(hits) > 1:
        logger.warning("%d rows match %s; using the first", len(hits), variant_id)
    row = df.loc[hits[0]]
    records = []
    group_order = []
    for group in pd.unique(sample_info["group"].astype(str)):
        chosen = _samples_per_biorep(sample_info, group)
        if group not in group_order:
            group_order.append(group)
        for _, info in chosen.iterrows():
            col = f"{info['sample_id']}_RPM"
            if col not in df.columns:
                continue
            records.append({"group": group, "bio_rep": info["bio_rep"], "rpm": float(row[col])})
    if not records:
        logger.warning("No RPM columns available for %s", variant_id)
        return False
    long_df = pd.DataFrame(records)
    summary = (
        long_df.groupby("group", sort=False)["rpm"]
        .agg(mean="mean", sd=lambda s: float(s.std(ddof=1)) if s.notna().sum() >= 2 else np.nan, n="count")
        .reindex(group_order)
        .reset_index()
    )
    fig, ax = plt.subplots(figsize=(7.5, 5))
    x = np.arange(len(summary))
    ax.bar(x, summary["mean"].to_numpy(dtype=float), color="#4575b4", width=0.65, alpha=0.92)
    for i, rec in summary.iterrows():
        if rec["n"] >= 2 and np.isfinite(rec["sd"]):
            ymin = max(0.0, float(rec["mean"]) - float(rec["sd"]))
            ymax = float(rec["mean"]) + float(rec["sd"])
            ax.errorbar(i, rec["mean"], yerr=[[rec["mean"] - ymin], [ymax - rec["mean"]]], fmt="none", ecolor="black", elinewidth=0.6, capsize=2)
    rng = np.random.default_rng(0)
    for i, group in enumerate(group_order):
        pts = long_df.loc[long_df["group"] == group, "rpm"].to_numpy(dtype=float)
        jitter = rng.uniform(-0.12, 0.12, size=len(pts))
        ax.scatter(np.full(len(pts), i) + jitter, pts, s=22, c="#2c2c2c", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(group_order, rotation=35, ha="right")
    ax.set_ylabel("RPM")
    ax.set_title(f"{variant_id} RPM (contamination)")
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    logger.info("Saved %s RPM bar plot to %s", variant_id, out_pdf)
    return True


def _integer_counts(raw: pd.DataFrame) -> pd.DataFrame:
    numeric = raw.apply(pd.to_numeric, errors="coerce").fillna(0)
    rounded = numeric.round()
    if ((numeric - rounded).abs() > 1e-6).any().any():
        logger.warning("Rounding non-integer counts before DESeq2")
    rounded = rounded.astype(int)
    if (rounded < 0).any().any():
        raise ValueError("Count matrix contains negative values")
    return rounded


def build_deseq_design(
    sample_info: pd.DataFrame,
    adjust_batch: bool = True,
    batch_col: str = "batch_seq",
) -> Tuple[str, List[str]]:
    """
    Choose the DESeq2 design for CapScreen samples.

    Always include ``group``. Never include ``bio_rep``: replicate numbers are
    labels within a group, not a shared blocking factor across groups. Include
    ``batch_seq`` only when ``adjust_batch`` is true, at least two batches are
    present, and batch is not fully confounded with group (each group exclusive
    to one batch and each batch exclusive to one group).

    Returns
    -------
    design :
        Formula string, e.g. ``~group`` or ``~batch_seq + group``.
    metadata_cols :
        Columns required in the colData for that design.
    """
    metadata_cols = ["group"]
    if not adjust_batch:
        return "~group", metadata_cols
    if batch_col not in sample_info.columns:
        logger.info("No %s column in sample_info; using design ~group", batch_col)
        return "~group", metadata_cols

    batches = sample_info[batch_col]
    n_batches = int(pd.Series(batches).nunique(dropna=True))
    if n_batches < 2:
        logger.info("Single sequencing batch detected; using design ~group")
        return "~group", metadata_cols

    groups = sample_info["group"].astype(str)
    batch_labels = batches.astype(str)
    n_batch_per_group = (
        pd.DataFrame({"group": groups, "batch": batch_labels})
        .groupby("group")["batch"]
        .nunique()
    )
    n_group_per_batch = (
        pd.DataFrame({"group": groups, "batch": batch_labels})
        .groupby("batch")["group"]
        .nunique()
    )
    if bool((n_batch_per_group == 1).all() and (n_group_per_batch == 1).all()):
        logger.warning(
            "batch_seq is fully confounded with group "
            "(each group is exclusive to one batch). "
            "Cannot adjust for batch; using design ~group."
        )
        return "~group", metadata_cols

    logger.info(
        "Multiple sequencing batches detected (%d); using design ~%s + group",
        n_batches,
        batch_col,
    )
    return f"~{batch_col} + group", [batch_col, "group"]


def _deseq_metadata(
    sample_info: pd.DataFrame,
    sample_ids: List[str],
    metadata_cols: List[str],
) -> pd.DataFrame:
    """Build DESeq2 colData with categorical factors for the chosen design."""
    missing = [col for col in metadata_cols if col not in sample_info.columns]
    if missing:
        raise ValueError(f"sample_info is missing design columns: {missing}")
    metadata = sample_info.loc[sample_ids, metadata_cols].copy()
    for col in metadata_cols:
        metadata[col] = metadata[col].astype(str)
    return metadata


def fit_blind_vst(
    counts: pd.DataFrame,
    sample_info: pd.DataFrame,
    n_cpus: int,
    adjust_batch: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build a DESeqDataSet and return blind VST counts.

    Blind VST fits an intercept-only dispersion trend (``use_design=False``).
    The stored design still follows :func:`build_deseq_design` so batch is
    recorded consistently when present.
    """
    sample_ids = [sid for sid in sample_info["sample_id"].astype(str) if sid in counts.columns]
    if len(sample_ids) < 2:
        raise ValueError("Need at least 2 samples with count columns for VST")
    if "ID" not in counts.columns:
        raise ValueError("Count table must contain ID for VST")
    count_mat = _integer_counts(counts[sample_ids])
    count_mat.index = counts["ID"].astype(str).to_numpy()
    if count_mat.index.duplicated().any():
        raise ValueError("Duplicate ID values while building the VST dataset")
    count_mat = count_mat.T
    count_mat.index = sample_ids
    design, metadata_cols = build_deseq_design(sample_info.loc[sample_ids], adjust_batch=adjust_batch)
    metadata = _deseq_metadata(sample_info, sample_ids, metadata_cols)
    dds = DeseqDataSet(
        counts=count_mat,
        metadata=metadata,
        design=design,
        refit_cooks=True,
        n_cpus=n_cpus,
        quiet=True,
    )
    dds.vst(use_design=False)
    vst = pd.DataFrame(dds.layers["vst_counts"], index=dds.obs_names, columns=dds.var_names).T
    vst.index.name = "ID"
    return vst, metadata


def drop_control_variants(vst: pd.DataFrame, pattern: str = "AAV") -> pd.DataFrame:
    """Drop control variants such as AAV2 from a variants-by-samples matrix."""
    if vst.empty:
        return vst
    mask = vst.index.astype(str).str.contains(pattern, case=False, na=False)
    n_drop = int(mask.sum())
    if n_drop:
        logger.info("Excluding %d control variant(s) matching %r from PCA", n_drop, pattern)
        return vst.loc[~mask].copy()
    return vst


def plot_vst_pearson_heatmap(vst: pd.DataFrame, out_pdf: Path) -> pd.DataFrame:
    """Sample-sample Pearson correlation of blind VST values.

    Color limits follow the observed correlation range so weak replicate
    correlations are not compressed into the same color as 0.7.
    """
    corr = vst.corr(method="pearson")
    finite = corr.to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise ValueError("No finite Pearson correlations to plot")
    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    if np.isclose(vmin, vmax):
        pad = 0.01 if vmax < 1.0 else 0.0
        vmin = max(-1.0, vmin - pad)
        vmax = min(1.0, vmax + pad)
    cmap = LinearSegmentedColormap.from_list(
        "vst_corr",
        ["#313695", "#4575b4", "#74add1", "#abd9e9", "#f7f7f7", "#fdae61", "#f46d43", "#d73027", "#a50026"],
    )
    n = corr.shape[0]
    fig_w = max(8.5, 0.35 * n + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_w * 0.9))
    sns.heatmap(
        corr,
        ax=ax,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        square=True,
        cbar_kws={"label": "Pearson r"},
    )
    ax.set_title("Sample Pearson correlation (blind VST)")
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    logger.info(
        "Saved VST Pearson heatmap to %s (color scale %.3f to %.3f)",
        out_pdf,
        vmin,
        vmax,
    )
    return corr


def _pca_frame(vst: pd.DataFrame, sample_info: pd.DataFrame, n_top: Optional[int], exclude_groups: List[str]):
    meta = sample_info.copy()
    meta["group"] = meta["group"].astype(str)
    keep = [sid for sid in vst.columns if sid in meta.index and meta.loc[sid, "group"] not in set(exclude_groups)]
    if len(keep) < 2:
        raise ValueError(
            f"Need at least 2 samples after excluding {exclude_groups}; remaining: {keep}"
        )
    mat = drop_control_variants(vst[keep].apply(pd.to_numeric, errors="coerce"))
    if mat.empty:
        raise ValueError("No variants left for PCA after excluding control variants")
    variance = mat.var(axis=1, ddof=1)
    ranked = variance[(variance > 0) & np.isfinite(variance)].sort_values(ascending=False)
    if ranked.empty:
        raise ValueError("No variants with positive finite VST variance")
    if n_top is None:
        selected = ranked.index
    else:
        selected = ranked.index[: int(n_top)]
    features = mat.loc[selected].T
    n_components = min(2, features.shape[0], features.shape[1])
    pca = PCA(n_components=n_components, svd_solver="full")
    coords = pca.fit_transform(features)
    frame = pd.DataFrame(coords, index=features.index, columns=[f"PC{i + 1}" for i in range(n_components)])
    if "PC2" not in frame.columns:
        frame["PC2"] = 0.0
    frame["group"] = meta.loc[frame.index, "group"].astype(str).to_numpy()
    frame["bio_rep"] = meta.loc[frame.index, "bio_rep"].astype(str).to_numpy()
    frame["sample"] = frame.index.astype(str)
    pct = pca.explained_variance_ratio_ * 100.0
    return frame, pct, len(selected)


def _draw_pca(
    frame: pd.DataFrame,
    pct: np.ndarray,
    n_variants: int,
    title: str,
    exclude_groups: List[str],
    fig: plt.Figure,
    ax: plt.Axes,
    point_size: float = 120,
    fontsize: float = 10,
) -> None:
    groups = list(dict.fromkeys(frame["group"]))
    markers = ["o", "^", "s", "P", "X", "*", "D", "v", "<", ">", "p", "h"]
    marker_map = {group: markers[i % len(markers)] for i, group in enumerate(groups)}
    sns.scatterplot(
        data=frame,
        x="PC1",
        y="PC2" if "PC2" in frame.columns else "PC1",
        hue="bio_rep",
        style="group",
        markers=marker_map,
        s=point_size,
        ax=ax,
    )
    pc1 = f"PC1 ({pct[0]:.0f}% variance)" if len(pct) else "PC1"
    pc2 = f"PC2 ({pct[1]:.0f}% variance)" if len(pct) > 1 else "PC2"
    ax.set_xlabel(pc1, fontsize=fontsize)
    ax.set_ylabel(pc2, fontsize=fontsize)
    excluded = ", ".join(exclude_groups) if exclude_groups else "none"
    subtitle = f"n variants = {n_variants} | excluded groups: {excluded}"
    # "outside" legend locs require a figure legend (matplotlib >= 3.6).
    handles, labels = ax.get_legend_handles_labels()
    if ax.get_legend() is not None:
        ax.get_legend().remove()
    fig.legend(
        handles,
        labels,
        loc="outside right upper",
        fontsize=max(fontsize - 2, 6),
        frameon=True,
    )
    _set_title_subtitle(ax, title, subtitle)


def plot_vst_pca(
    vst: pd.DataFrame,
    sample_info: pd.DataFrame,
    qc_dir: Path,
    n_top: int,
    exclude_groups: List[str],
    point_size: float = 120,
    fontsize: float = 10,
    figsize: Optional[Tuple[float, float]] = None,
) -> None:
    """Write the two notebook PCA figures: all positive-variance variants, and the top variable set.

    Input samples are dropped via ``exclude_groups``. Control variants matching
    AAV (for example AAV2) are removed so contamination does not drive the PCs.
    """
    all_frame, all_pct, n_all = _pca_frame(vst, sample_info, None, exclude_groups)
    top_frame, top_pct, n_top_used = _pca_frame(vst, sample_info, n_top, exclude_groups)
    fig_w, fig_h = figsize if figsize is not None else (9.0, 7.0)
    plots = [
        (all_frame, all_pct, n_all, "PCA on VST (all variants with positive variance)", qc_dir / "pca_all_features.png"),
        (top_frame, top_pct, n_top_used, f"PCA on VST (top {n_top} variable variants)", qc_dir / "pca_top500_features.png"),
    ]
    with PdfPages(qc_dir / "pca_vst.pdf") as pdf:
        for frame, pct, n_variants, title, png_path in plots:
            fig, ax = plt.subplots(figsize=(fig_w, fig_h), layout="constrained")
            _draw_pca(
                frame, pct, n_variants, title, exclude_groups, fig, ax,
                point_size=point_size, fontsize=fontsize,
            )
            fig.savefig(png_path, dpi=300)
            pdf.savefig(fig)
            plt.close(fig)
            logger.info("Saved %s", png_path)
    logger.info("Saved two-page VST PCA to %s", qc_dir / "pca_vst.pdf")


def _samples_for_group(
    sample_info: pd.DataFrame,
    sample_ids: List[str],
    group: str,
) -> List[str]:
    return [
        sid
        for sid in sample_ids
        if str(sample_info.loc[sid, "group"]) == str(group)
    ]


def _finalize_contrast_results(
    results_df: pd.DataFrame,
    raw_counts: pd.DataFrame,
    all_variant_ids: pd.Index,
    y_samples: List[str],
    x_samples: List[str],
    y_group: str,
    x_group: str,
    design: str,
) -> pd.DataFrame:
    """
    Reindex DESeq2 results to every filtered variant and blank unsupported tests.

    Variants with zero counts in both contrast groups cannot support a fold
    change. Keep them in the table with missing test statistics (and baseMean 0),
    matching the R notebook behavior.
    """
    res = results_df.copy()
    res.index = res.index.astype(str)
    res = res.reindex(all_variant_ids.astype(str))
    res.insert(0, "ID", all_variant_ids.astype(str).to_numpy())

    y_sum = raw_counts.loc[all_variant_ids, y_samples].sum(axis=1)
    x_sum = raw_counts.loc[all_variant_ids, x_samples].sum(axis=1)
    zero_both = (y_sum.to_numpy() == 0) & (x_sum.to_numpy() == 0)
    n_zero = int(zero_both.sum())
    if n_zero:
        for col in ("log2FoldChange", "lfcSE", "stat", "pvalue", "padj"):
            if col in res.columns:
                res.loc[zero_both, col] = np.nan
        if "baseMean" in res.columns:
            res.loc[zero_both, "baseMean"] = 0.0
        logger.info(
            "%s vs %s: %d variant(s) have zero counts in both groups; "
            "wrote missing DE statistics",
            y_group,
            x_group,
            n_zero,
        )

    res["y_group"] = y_group
    res["x_group"] = x_group
    res["comparison"] = f"{y_group}_vs_{x_group}"
    res["design"] = design
    return res


def run_deseq2_comparisons(
    counts: pd.DataFrame,
    sample_info: pd.DataFrame,
    comparisons: pd.DataFrame,
    out_dir: Path,
    n_cpus: int,
    min_replicates: int = 2,
    adjust_batch: bool = True,
    write_excel: bool = False,
) -> int:
    """
    Fit each contrast pairwise on only the samples in its two groups.

    Size factors and dispersions are estimated per contrast from that subset
    alone, so results do not depend on other groups in the experiment (e.g.
    low-variance Input). Design comes from :func:`build_deseq_design` on the
    subset (``~group`` or ``~batch_seq + group``). ``bio_rep`` is never in the
    design.

    Variants with zero counts in both groups of a contrast are kept in the
    output table with missing test statistics.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    all_sample_ids = [
        sid for sid in sample_info["sample_id"].astype(str) if sid in counts.columns
    ]
    if len(all_sample_ids) < 2:
        raise ValueError("Need at least 2 samples with count columns for DESeq2")
    if "ID" not in counts.columns:
        raise ValueError("Count table must contain ID for DESeq2")

    raw_all = _integer_counts(counts[all_sample_ids])
    raw_all.index = counts["ID"].astype(str)
    if raw_all.index.duplicated().any():
        raise ValueError("Duplicate variant IDs while preparing the DESeq2 dataset")
    all_variant_ids = raw_all.index

    group_counts = sample_info.loc[all_sample_ids, "group"].astype(str).value_counts()
    n_ok = 0
    for _, cmp in comparisons.iterrows():
        y_group = str(cmp["y_group"])
        x_group = str(cmp["x_group"])
        name = f"{y_group}_vs_{x_group}"
        n_y = int(group_counts.get(y_group, 0))
        n_x = int(group_counts.get(x_group, 0))
        if n_y < min_replicates or n_x < min_replicates:
            logger.warning(
                "Skipping %s: need >= %d sample(s) per group; counts: %s=%d, %s=%d",
                name,
                min_replicates,
                y_group,
                n_y,
                x_group,
                n_x,
            )
            continue

        y_samples = _samples_for_group(sample_info, all_sample_ids, y_group)
        x_samples = _samples_for_group(sample_info, all_sample_ids, x_group)
        sample_ids = y_samples + x_samples
        if len(sample_ids) < 2:
            logger.warning("Skipping %s: fewer than 2 samples after group selection", name)
            continue

        subset_info = sample_info.loc[sample_ids]
        count_mat = raw_all[sample_ids]
        keep_genes = count_mat.sum(axis=1) > 0
        n_drop = int((~keep_genes).sum())
        if n_drop:
            logger.info(
                "%s: omitting %d variant(s) with zero counts in the contrasted "
                "samples from the fit; they remain in the table with missing statistics",
                name,
                n_drop,
            )
        count_mat = count_mat.loc[keep_genes]
        if count_mat.empty:
            logger.warning("Skipping %s: no expressed variants in the contrasted samples", name)
            continue

        design, metadata_cols = build_deseq_design(subset_info, adjust_batch=adjust_batch)
        metadata = _deseq_metadata(subset_info, sample_ids, metadata_cols)
        design_counts = count_mat.T

        logger.info(
            "Fitting pairwise DESeq2 for %s on %d samples x %d variants (design=%s)",
            name,
            design_counts.shape[0],
            design_counts.shape[1],
            design,
        )
        try:
            dds = DeseqDataSet(
                counts=design_counts,
                metadata=metadata,
                design=design,
                refit_cooks=True,
                n_cpus=n_cpus,
                quiet=True,
            )
            dds.deseq2()
            stats = DeseqStats(dds, contrast=["group", y_group, x_group], quiet=True)
            stats.summary()
            res = _finalize_contrast_results(
                stats.results_df,
                raw_all,
                all_variant_ids,
                y_samples,
                x_samples,
                y_group,
                x_group,
                design,
            )
            out_path = out_dir / de_result_filename(y_group, x_group)
            write_table(res, out_path, write_excel=write_excel)
            n_sig = int((res["padj"] < 0.05).fillna(False).sum()) if "padj" in res.columns else 0
            logger.info(
                "Wrote %s (pairwise, design=%s, %d variants, %d with padj < 0.05)",
                out_path,
                design,
                len(res),
                n_sig,
            )
            n_ok += 1
        except Exception as exc:
            logger.warning("DESeq2 contrast failed for %s: %s", name, exc)
    return n_ok


def _cfg_int(value, default: int) -> int:
    if value is None:
        return default
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


_EXCEL_TEXT_COLS = frozenset({
    "ID", "ID_WLG", "peptide", "barcode", "variable_seq",
    "y_group", "x_group", "comparison", "design", "arm", "tag",
    "target", "backgrounds", "input_group", "title", "subtitle",
})


def write_table(df: pd.DataFrame, path: Path, *, write_excel: bool = False) -> None:
    """
    Write a results table as TSV/CSV and optionally as Excel.

    TSV/CSV remains the canonical machine-readable output. When
    ``write_excel`` is true, also write ``path.with_suffix('.xlsx')``.
    ID-like columns are forced to text in Excel to avoid mangling.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".tsv":
        df.to_csv(path, sep="\t", index=False)
    else:
        df.to_csv(path, index=False)

    if not write_excel:
        return

    xlsx_path = path.with_suffix(".xlsx")
    try:
        out = df.copy()
        for col in out.columns:
            if col in _EXCEL_TEXT_COLS or str(col).endswith("_id"):
                out[col] = out[col].astype(str).where(out[col].notna(), "")
        out.to_excel(xlsx_path, index=False, engine="openpyxl")
        logger.info("Wrote Excel companion %s", xlsx_path)
    except ImportError:
        logger.warning(
            "write_excel is enabled but openpyxl is not installed; "
            "skipped %s (pip/conda install openpyxl in the capscreen-stat env)",
            xlsx_path,
        )
    except Exception as exc:
        logger.warning("Failed to write Excel companion %s: %s", xlsx_path, exc)


def _rpm_cols_for_group(
    sample_info: pd.DataFrame,
    group: str,
    df: pd.DataFrame,
) -> List[str]:
    """Return ``{sample_id}_RPM`` columns for a sample_info group."""
    ids = (
        sample_info.loc[sample_info["group"].astype(str) == str(group), "sample_id"]
        .astype(str)
        .tolist()
    )
    cols = [f"{sid}_RPM" for sid in ids if f"{sid}_RPM" in df.columns]
    if not cols:
        raise ValueError(f"No RPM columns found for group '{group}'")
    return cols


def _row_median(mat: np.ndarray) -> np.ndarray:
    return np.nanmedian(mat, axis=1)


def _er_eh_ratio(num: np.ndarray, den: np.ndarray, pseudo: float = 1.0) -> np.ndarray:
    """Lin et al. abundance ratio with a 1 RPM pseudocount on both arms."""
    return (num + pseudo) / (den + pseudo)


def _variant_display_label(variant_id: str) -> str:
    text = str(variant_id)
    if text.startswith("WID-WID"):
        return text[len("WID-WID"):]
    if text.startswith("WID-"):
        return text[len("WID-"):]
    return text


def _model_key(group: str) -> str:
    """Leading cohort/model token (OA, OPA, FASMA, …) used to scope RPM bars."""
    text = str(group)
    for marker in ("-WT-", "-KO-"):
        if marker in text:
            return text.split(marker, 1)[0]
    if "-" in text:
        return text.split("-", 1)[0]
    return text


def _bar_groups_for_arm(
    sample_info: pd.DataFrame,
    target: str,
    backgrounds: List[str],
    input_group: str,
) -> List[str]:
    """
    Groups for RPM bars, matching the R arm helpers.

    Keep Input, the arm target/backgrounds, and other Sample_info groups that
    share the target's model/cohort key (e.g. all ``OPA-*`` for an OPA arm).
    """
    ordered = _sample_info_groups_in_order(sample_info)
    key = _model_key(target)
    input_names = {str(input_group), "Input"}
    required = {str(target), *[str(b) for b in backgrounds]}
    keep = []
    for group in ordered:
        g = str(group)
        if g in input_names or g in required or _model_key(g) == key:
            keep.append(g)
    return keep


def _set_title_subtitle(ax: plt.Axes, title: str, subtitle: str) -> None:
    """Match R labs(title, subtitle): bold ~12pt title, ~8.5pt subtitle, left."""
    title_s = str(title).strip()
    sub_lines = [ln.strip() for ln in str(subtitle).split("\n") if ln.strip()]
    # Matplotlib has no native subtitle; pack under one left title. Constrained
    # layout then reserves a tight band so the legend does not collide with it.
    if sub_lines:
        ax.set_title(
            title_s + "\n" + "\n".join(sub_lines),
            loc="left",
            fontsize=10,
            fontweight="bold",
            linespacing=1.35,
            pad=6,
            color="0.1",
        )
    else:
        ax.set_title(
            title_s,
            loc="left",
            fontsize=12,
            fontweight="bold",
            pad=6,
            color="0.1",
        )


def _apply_log10_ticks(ax: plt.Axes, *, which: str = "both") -> None:
    """Major 10^n tick labels; suppress noisy minor labels."""
    formatter = LogFormatterMathtext(base=10)
    locator = LogLocator(base=10)
    if which in ("both", "x"):
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_minor_locator(NullLocator())
    if which in ("both", "y"):
        ax.yaxis.set_major_locator(LogLocator(base=10))
        ax.yaxis.set_major_formatter(formatter)
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs="auto"))


def _sample_info_groups_in_order(
    sample_info: pd.DataFrame,
    group_col: str = "group",
) -> List[str]:
    """Unique group labels in first-seen Sample_info order."""
    return list(dict.fromkeys(sample_info[group_col].astype(str).tolist()))


def _validate_enrichment_metadata_consistency(
    sample_info: pd.DataFrame,
    group_col: str = "group",
) -> pd.DataFrame:
    """Normalize enrichment columns and require consistent values within a group."""
    info = sample_info.copy()
    missing = [c for c in ("enrichment_role", "enrichment_arm") if c not in info.columns]
    if missing:
        raise ValueError(
            "Sample_info must define DESeq2 and enrichment/enhancement contrasts via "
            f"columns {missing}. Required: enrichment_role "
            "(input|target|background) and enrichment_arm. See README."
        )
    info["_role"] = info["enrichment_role"].fillna("").astype(str).str.strip().str.lower()
    info["_arm"] = info["enrichment_arm"].fillna("").astype(str).str.strip()
    if not (info["_role"] != "").any():
        raise ValueError(
            "Sample_info enrichment_role is empty for all samples; assign at least "
            "one input, one target, and one background group"
        )

    allowed = {"", "input", "target", "background"}
    bad = sorted(set(info["_role"]) - allowed)
    if bad:
        raise ValueError(
            "enrichment_role must be input, target, background, or blank; got: "
            + ", ".join(bad)
        )

    for group, sub in info.groupby(group_col, sort=False):
        roles = sorted(set(sub["_role"]))
        arms = sorted(set(sub["_arm"]))
        if len(roles) > 1:
            raise ValueError(
                f"Group '{group}' has inconsistent enrichment_role values: {roles}"
            )
        if len(arms) > 1:
            raise ValueError(
                f"Group '{group}' has inconsistent enrichment_arm values: {arms}"
            )
        role = roles[0]
        arm = arms[0]
        if role in {"target", "background"} and not arm:
            raise ValueError(
                f"Group '{group}' has enrichment_role={role} but empty enrichment_arm"
            )
        if role == "input" and arm:
            logger.warning(
                "Input group '%s' has enrichment_arm=%r; arm is ignored for input",
                group,
                arm,
            )
    return info


def _arms_from_enrichment_metadata(
    sample_info: pd.DataFrame,
    group_col: str = "group",
) -> pd.DataFrame:
    """
    Build Er/Eh arms from Sample_info enrichment_role / enrichment_arm.

    enrichment_role: input | target | background
    enrichment_arm: free-form tag linking one target group to its background
    groups. A target may list several arms separated by ``;``.
    """
    info = _validate_enrichment_metadata_consistency(sample_info, group_col=group_col)
    input_groups = list(
        dict.fromkeys(
            info.loc[info["_role"] == "input", group_col].astype(str).tolist()
        )
    )
    if not input_groups:
        raise ValueError(
            "Sample_info enrichment_role is set but no group has enrichment_role=input"
        )
    if len(input_groups) > 1:
        raise ValueError(
            "Sample_info defines multiple input groups via enrichment_role=input: "
            + ", ".join(input_groups)
        )
    resolved_input = input_groups[0]

    def _split_arms(value: str) -> List[str]:
        return [a.strip() for a in str(value).split(";") if a.strip()]

    tb = info.loc[info["_role"].isin(["target", "background"])].copy()
    arm_names: List[str] = []
    for raw in tb["_arm"].tolist():
        for arm_tag in _split_arms(raw):
            if arm_tag not in arm_names:
                arm_names.append(arm_tag)
    if not arm_names:
        raise ValueError(
            "Sample_info enrichment_role is set but no enrichment_arm values were found "
            "for target/background groups"
        )

    arms: List[Dict[str, Any]] = []
    for arm_tag in arm_names:
        members = tb.loc[tb["_arm"].map(lambda v: arm_tag in _split_arms(v))]
        targets = list(
            dict.fromkeys(members.loc[members["_role"] == "target", group_col].astype(str))
        )
        backgrounds = list(
            dict.fromkeys(
                members.loc[members["_role"] == "background", group_col].astype(str)
            )
        )
        if len(targets) == 0:
            raise ValueError(f"enrichment_arm '{arm_tag}' has no target group")
        if len(targets) > 1:
            raise ValueError(
                f"enrichment_arm '{arm_tag}' has multiple target groups: {', '.join(targets)}"
            )
        if not backgrounds:
            raise ValueError(f"enrichment_arm '{arm_tag}' has no background group")
        target = targets[0]
        file_tag = sanitize_group_label(arm_tag)
        arms.append(
            {
                "tag": file_tag,
                "target": target,
                "backgrounds": ";".join(backgrounds),
                "input_group": resolved_input,
                "title": f"Enrichment/enhancement arm {arm_tag}",
                "subtitle": (
                    f"Er = {target} / {resolved_input}; "
                    f"Eh = {target} / {' + '.join(backgrounds)}"
                ),
                "mode": "sample_info",
            }
        )
    out = pd.DataFrame(arms)
    logger.info(
        "Built %d enrichment/enhancement arm(s) from Sample_info enrichment_role/arm "
        "(input=%s)",
        len(out),
        resolved_input,
    )
    return out


def build_enrichment_enhancement_arms(
    sample_info: pd.DataFrame,
    group_col: str = "group",
) -> pd.DataFrame:
    """
    Build Lin-style Er/Eh arms only from Sample_info roles/arms.

    Requires ``enrichment_role`` / ``enrichment_arm``. Target versus background,
    with the Sample_info input group as the Er denominator.
    """
    return _arms_from_enrichment_metadata(sample_info, group_col=group_col)


def score_selection(
    df: pd.DataFrame,
    sample_info: pd.DataFrame,
    target: str,
    backgrounds: List[str],
    input_group: str = "Input",
    pseudocount: float = 1.0,
) -> pd.DataFrame:
    """
    Score one selection arm against its input pool and background arm(s).

    Er = A_i(target) / A_i(input); Eh = A_i(target) / A_i(backgrounds).
    A_i is RPM. Both replicate-median and replicate-mean summaries are reported,
    plus per-target-replicate Er/Eh against the arm medians.
    """
    tgt_cols = _rpm_cols_for_group(sample_info, target, df)
    bg_cols: List[str] = []
    for bg in backgrounds:
        bg_cols.extend(_rpm_cols_for_group(sample_info, bg, df))
    inp_cols = _rpm_cols_for_group(sample_info, input_group, df)

    ids = df["ID"].astype(str).to_numpy()
    peptide = (
        df["peptide"].astype(str).to_numpy()
        if "peptide" in df.columns
        else np.array([""] * len(df), dtype=object)
    )
    tgt = df[tgt_cols].to_numpy(dtype=float)
    bg = df[bg_cols].to_numpy(dtype=float)
    inp = df[inp_cols].to_numpy(dtype=float)

    tgt_med = _row_median(tgt)
    bg_med = _row_median(bg)
    inp_med = _row_median(inp)
    tgt_mean = np.nanmean(tgt, axis=1)
    bg_mean = np.nanmean(bg, axis=1)
    inp_mean = np.nanmean(inp, axis=1)

    er_reps = _er_eh_ratio(tgt, inp_med[:, None], pseudocount)
    eh_reps = _er_eh_ratio(tgt, bg_med[:, None], pseudocount)

    out = pd.DataFrame(
        {
            "ID": ids,
            "peptide": peptide,
            "is_control": ~pd.Series(ids).str.startswith("WID-", na=False).to_numpy(),
            "target_RPM_median": tgt_med,
            "background_RPM_median": bg_med,
            "input_RPM_median": inp_med,
            "target_RPM_mean": tgt_mean,
            "background_RPM_mean": bg_mean,
            "input_RPM_mean": inp_mean,
            "Er_median": _er_eh_ratio(tgt_med, inp_med, pseudocount),
            "Eh_median": _er_eh_ratio(tgt_med, bg_med, pseudocount),
            "Er_mean": _er_eh_ratio(tgt_mean, inp_mean, pseudocount),
            "Eh_mean": _er_eh_ratio(tgt_mean, bg_mean, pseudocount),
        }
    )
    for i in range(er_reps.shape[1]):
        out[f"Er_rep{i + 1}"] = er_reps[:, i]
        out[f"Eh_rep{i + 1}"] = eh_reps[:, i]
    return out


def call_top_er_eh(
    scores: pd.DataFrame,
    estimator: str = "median",
    er_cut: float = 3.0,
    eh_cut: float = 3.0,
    min_rpm: float = 5.0,
    n_top: Optional[int] = 10,
) -> pd.DataFrame:
    """Call top performers inside the Er/Eh performance box."""
    er_rep_cols = [c for c in scores.columns if re.fullmatch(r"Er_rep\d+", c)]
    eh_rep_cols = [c for c in scores.columns if re.fullmatch(r"Eh_rep\d+", c)]
    if not er_rep_cols or not eh_rep_cols:
        raise ValueError("score table is missing per-replicate Er/Eh columns")
    n_reps = len(er_rep_cols)

    out = scores.copy()
    out["Er"] = out[f"Er_{estimator}"]
    out["Eh"] = out[f"Eh_{estimator}"]
    out["target_RPM"] = out[f"target_RPM_{estimator}"]
    out["n_reps_pass"] = (
        (out[er_rep_cols].to_numpy(dtype=float) > er_cut)
        & (out[eh_rep_cols].to_numpy(dtype=float) > eh_cut)
    ).sum(axis=1)
    out["confidence"] = out["n_reps_pass"] / float(n_reps)
    out["estimator"] = estimator

    hits = out.loc[
        (~out["is_control"])
        & (out["Er"] > er_cut)
        & (out["Eh"] > eh_cut)
        & (out["target_RPM"] >= min_rpm)
        & (out["n_reps_pass"] >= int(np.ceil(n_reps / 2.0)))
    ].copy()
    score = np.log2(hits["Er"].to_numpy(dtype=float)) + np.log2(hits["Eh"].to_numpy(dtype=float))
    hits = hits.assign(_score=score).sort_values(
        by=["confidence", "_score"], ascending=[False, False]
    ).drop(columns=["_score"])
    hits["rank"] = np.arange(1, len(hits) + 1)
    if n_top is not None:
        hits = hits.head(int(n_top)).copy()
    return hits.reset_index(drop=True)


def plot_er_eh_scatter(
    scores: pd.DataFrame,
    top: pd.DataFrame,
    estimator: str,
    out_prefix: Path,
    title: str,
    subtitle: str,
    er_cut: float = 3.0,
    eh_cut: float = 3.0,
    dpi: int = 300,
) -> None:
    """Paper-style log-log enrichment versus enhancement scatter."""
    d = scores.copy()
    d["Er"] = d[f"Er_{estimator}"]
    d["Eh"] = d[f"Eh_{estimator}"]
    bg = d.loc[~d["is_control"]]
    ctrl = d.loc[d["is_control"]]

    vals = np.concatenate([d["Er"].to_numpy(dtype=float), d["Eh"].to_numpy(dtype=float)])
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size == 0:
        logger.warning("No finite Er/Eh values for %s", out_prefix.name)
        return
    lo = float(vals.min())
    hi = float(vals.max()) * 1.6
    lims = (max(lo, 1e-3), hi)

    n_top = 0 if top is None or top.empty else len(top)
    # Match R enrichment_enhancement.R ggsave(width=9.5, height=6.5).
    fig, ax = plt.subplots(figsize=(9.5, 6.5), layout="constrained")
    ax.scatter(
        bg["Er"], bg["Eh"],
        s=6, c="0.72", alpha=0.35, linewidths=0, zorder=1, rasterized=True,
    )
    ax.plot(lims, lims, ls="--", color="0.45", lw=1, zorder=0)
    if er_cut < lims[1] and eh_cut < lims[1]:
        ax.add_patch(
            plt.Rectangle(
                (er_cut, eh_cut),
                lims[1] - er_cut,
                lims[1] - eh_cut,
                fill=False,
                edgecolor="goldenrod",
                linestyle=":",
                linewidth=1.0,
                zorder=2,
            )
        )
    if not ctrl.empty:
        ax.scatter(ctrl["Er"], ctrl["Eh"], s=36, c="0.15", zorder=3)
        for _, row in ctrl.iterrows():
            ax.annotate(
                str(row["ID"]),
                (row["Er"], row["Eh"]),
                textcoords="offset points",
                xytext=(4, 4),
                fontsize=7,
                fontweight="bold",
                color="0.15",
            )

    if n_top:
        hl = d.set_index("ID").loc[top["ID"]].reset_index()
        cmap = plt.get_cmap("Dark2")
        colors = [cmap(i % 8) for i in range(len(hl))]
        for i, row in hl.iterrows():
            peptide = row["peptide"] if "peptide" in hl.columns and pd.notna(row["peptide"]) else ""
            label = f"{row['ID']}\n{peptide}".strip()
            ax.scatter(
                [row["Er"]], [row["Eh"]],
                s=48, c=[colors[i]], zorder=4, label=label,
            )
        # R default legend is outside-right; figure legend required for "outside".
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            title=f"Top {len(hl)} reproducible\nvariants",
            fontsize=6.5,
            title_fontsize=7,
            loc="outside right upper",
            frameon=True,
            fancybox=False,
            edgecolor="0.75",
            borderpad=0.4,
            labelspacing=0.55,
            handletextpad=0.5,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lims)
    ax.set_ylim(*lims)
    ax.set_aspect("equal", adjustable="box")
    _apply_log10_ticks(ax)
    ax.set_xlabel("Enrichment (Er)")
    ax.set_ylabel("Enhancement (Eh)")
    ax.grid(True, which="major", linestyle="-", alpha=0.25)
    ax.grid(False, which="minor")
    _set_title_subtitle(ax, title, subtitle)

    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_prefix}.png", dpi=dpi)
    fig.savefig(f"{out_prefix}.pdf")
    plt.close(fig)


def plot_top_variant_rpm_bars(
    df: pd.DataFrame,
    sample_info: pd.DataFrame,
    ids: List[str],
    groups: List[str],
    out_prefix: Path,
    title: str,
    subtitle: str,
    dpi: int = 300,
) -> Optional[pd.DataFrame]:
    """
    Grouped RPM bar plot for ordered top variants.

    Bars are replicate medians; whiskers span the replicate range; individual
    replicates are overlaid. Linear and log10 y-axis versions are written.
    """
    if not ids:
        return None
    present_groups = []
    group_cols: Dict[str, List[str]] = {}
    for group in groups:
        try:
            cols = _rpm_cols_for_group(sample_info, group, df)
        except ValueError:
            continue
        present_groups.append(group)
        group_cols[group] = cols
    if not present_groups:
        logger.warning("No groups available for RPM bars (%s)", out_prefix.name)
        return None

    id_order = [str(v) for v in ids]
    sub = df.loc[df["ID"].astype(str).isin(id_order)].copy()
    if sub.empty:
        return None
    sub["_ord"] = sub["ID"].astype(str).map({v: i for i, v in enumerate(id_order)})
    sub = sub.sort_values("_ord")

    long_rows = []
    for group in present_groups:
        cols = group_cols[group]
        mat = sub[cols].to_numpy(dtype=float)
        for j, col in enumerate(cols):
            for i, vid in enumerate(sub["ID"].astype(str)):
                long_rows.append(
                    {
                        "ID": vid,
                        "variant": _variant_display_label(vid),
                        "group": group,
                        "replicate": j + 1,
                        "RPM": float(mat[i, j]),
                    }
                )
    long = pd.DataFrame(long_rows)
    long["variant"] = pd.Categorical(
        long["variant"],
        categories=[_variant_display_label(v) for v in id_order],
        ordered=True,
    )
    long["group"] = pd.Categorical(long["group"], categories=present_groups, ordered=True)

    stats = (
        long.groupby(["variant", "group"], observed=True)
        .agg(
            median_RPM=("RPM", "median"),
            mean_RPM=("RPM", "mean"),
            sd_RPM=("RPM", "std"),
            min_RPM=("RPM", "min"),
            max_RPM=("RPM", "max"),
            n=("RPM", "size"),
        )
        .reset_index()
    )

    n_groups = len(present_groups)
    cmap = plt.get_cmap("Dark2")
    colors = {g: cmap(i % 8) for i, g in enumerate(present_groups[:-1])}
    if present_groups:
        colors[present_groups[-1]] = "0.6"
        if present_groups[-1] != "Input" and "Input" in present_groups:
            colors["Input"] = "0.6"
            for i, g in enumerate([x for x in present_groups if x != "Input"]):
                colors[g] = cmap(i % 8)

    pos_rpm = long.loc[long["RPM"] > 0, "RPM"]
    floor_rpm = max(float(pos_rpm.min()) if len(pos_rpm) else 1e-3, 1e-3)
    variants = list(stats["variant"].cat.categories)
    x = np.arange(len(variants))
    width = 0.8 / max(n_groups, 1)

    def _draw(log_scale: bool) -> plt.Figure:
        # Match R enrichment_enhancement.R ggsave(width=11, height=5.5).
        ncol = min(n_groups, 5)
        fig, ax = plt.subplots(figsize=(11.0, 5.5), layout="constrained")
        for gi, group in enumerate(present_groups):
            gstat = stats.loc[stats["group"] == group].set_index("variant").reindex(variants)
            med = gstat["median_RPM"].to_numpy(dtype=float)
            ymin = gstat["min_RPM"].to_numpy(dtype=float)
            ymax = gstat["max_RPM"].to_numpy(dtype=float)
            # R uses pmax(value, floor_rpm) on log plots (bars grow from axis floor).
            if log_scale:
                med = np.where(np.isfinite(med) & (med > 0), np.maximum(med, floor_rpm), floor_rpm)
                ymin = np.where(np.isfinite(ymin) & (ymin > 0), np.maximum(ymin, floor_rpm), floor_rpm)
                ymax = np.where(np.isfinite(ymax) & (ymax > 0), np.maximum(ymax, floor_rpm), floor_rpm)
            else:
                med = np.nan_to_num(med, nan=0.0)
                ymin = np.nan_to_num(ymin, nan=0.0)
                ymax = np.nan_to_num(ymax, nan=0.0)
            xpos = x - 0.4 + width / 2 + gi * width
            ax.bar(
                xpos, med, width=width * 0.95,
                color=colors.get(group, "0.5"),
                edgecolor="0.25", linewidth=0.3,
                label=group, zorder=2,
            )
            yerr_lo = np.maximum(med - ymin, 0.0)
            yerr_hi = np.maximum(ymax - med, 0.0)
            ax.errorbar(
                xpos, med, yerr=[yerr_lo, yerr_hi],
                fmt="none", ecolor="0.25", elinewidth=0.6,
                capsize=2.0, capthick=0.6, zorder=3,
            )
            glong = long.loc[long["group"] == group]
            for vi, variant in enumerate(variants):
                pts = glong.loc[glong["variant"] == variant, "RPM"].to_numpy(dtype=float)
                if log_scale:
                    pts = np.where(np.isfinite(pts) & (pts > 0), np.maximum(pts, floor_rpm), floor_rpm)
                jitter = (np.random.default_rng(vi + gi * 17).random(len(pts)) - 0.5) * width * 0.35
                ax.scatter(
                    np.full(len(pts), xpos[vi]) + jitter,
                    pts,
                    s=12, facecolors="white", edgecolors="0.2",
                    linewidths=0.3, alpha=0.85, zorder=4,
                )

        ax.set_xticks(x)
        ax.set_xticklabels(variants, rotation=45, ha="right", fontsize=7.5)
        ax.set_xlabel("")
        ax.grid(False, axis="x")
        ax.grid(True, axis="y", which="major", alpha=0.25)
        ax.grid(False, which="minor")

        clip_note = False
        if log_scale:
            ax.set_yscale("log")
            _apply_log10_ticks(ax, which="y")
            ax.set_ylabel("RPM (replicate median, log10 scale)")
        else:
            ymax_plot = float(np.nanmax(stats["median_RPM"].to_numpy(dtype=float))) * 1.6
            if not np.isfinite(ymax_plot) or ymax_plot <= 0:
                ymax_plot = 1.0
            ax.set_ylim(0, ymax_plot)
            ax.set_ylabel("RPM (replicate median)")
            clip_note = bool((long["RPM"] > ymax_plot).any())

        _set_title_subtitle(ax, title, subtitle)
        # R theme: legend.position = "bottom" (enrichment_enhancement.R).
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="outside lower center",
            ncol=ncol,
            frameon=False,
            fontsize=8,
            columnspacing=1.0,
            handletextpad=0.4,
        )
        if clip_note:
            ax.set_xlabel(
                "Linear y axis clipped to the median range; some replicate values extend beyond it",
                fontsize=7.5, color="0.35", labelpad=2,
            )
        return fig

    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    stats.to_csv(f"{out_prefix}_data.csv", index=False)
    long.to_csv(f"{out_prefix}_replicates.csv", index=False)
    for name, is_log in (("linear", False), ("log", True)):
        fig = _draw(is_log)
        fig.savefig(f"{out_prefix}_{name}.png", dpi=dpi)
        fig.savefig(f"{out_prefix}_{name}.pdf")
        plt.close(fig)
    return stats


def run_enrichment_enhancement_arm(
    df: pd.DataFrame,
    sample_info: pd.DataFrame,
    arm: pd.Series,
    out_dir: Path,
    input_group: str = "Input",
    er_cut: float = 3.0,
    eh_cut: float = 3.0,
    min_rpm: float = 5.0,
    n_top: int = 10,
    pseudocount: float = 1.0,
    dpi: int = 300,
    write_excel: bool = False,
) -> Dict[str, Any]:
    """Score, call, and plot one enrichment/enhancement arm."""
    backgrounds = [c for c in str(arm["backgrounds"]).split(";") if c]
    tag = str(arm["tag"])
    target = str(arm["target"])
    title = str(arm["title"])
    subtitle = str(arm["subtitle"])
    arm_input = str(arm["input_group"]) if "input_group" in arm.index and pd.notna(arm["input_group"]) else input_group
    bar_groups = _bar_groups_for_arm(sample_info, target, backgrounds, arm_input)

    scores = score_selection(
        df, sample_info, target, backgrounds,
        input_group=arm_input, pseudocount=pseudocount,
    )
    hits = call_top_er_eh(
        scores, "median", er_cut=er_cut, eh_cut=eh_cut, min_rpm=min_rpm, n_top=None,
    )
    hits_mean = call_top_er_eh(
        scores, "mean", er_cut=er_cut, eh_cut=eh_cut, min_rpm=min_rpm, n_top=None,
    )
    top = hits.head(int(n_top)).copy()
    top_mean = hits_mean.head(int(n_top)).copy()

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_table(
        scores.sort_values("Er_median", ascending=False),
        out_dir / f"Er_Eh_{tag}_all_variants.csv",
        write_excel=write_excel,
    )
    write_table(hits, out_dir / f"Er_Eh_{tag}_hits_median.csv", write_excel=write_excel)
    write_table(hits_mean, out_dir / f"Er_Eh_{tag}_hits_mean.csv", write_excel=write_excel)
    write_table(top, out_dir / f"Er_Eh_{tag}_top_hits_median.csv", write_excel=write_excel)
    write_table(top_mean, out_dir / f"Er_Eh_{tag}_top_hits_mean.csv", write_excel=write_excel)

    plot_er_eh_scatter(
        scores, top, "median",
        out_dir / f"Er_Eh_{tag}_median",
        title,
        (
            f"{subtitle}\nA_i = replicate median RPM; highlighted variants clear "
            f"Er > {er_cut} and Eh > {eh_cut} in a majority of replicates"
        ),
        er_cut=er_cut, eh_cut=eh_cut, dpi=dpi,
    )
    plot_er_eh_scatter(
        scores, top, "mean",
        out_dir / f"Er_Eh_{tag}_mean",
        f"{title} (mean estimator)",
        (
            f"{subtitle}\nA_i = replicate mean RPM, as in the paper; same "
            "highlighted variants as the median figure"
        ),
        er_cut=er_cut, eh_cut=eh_cut, dpi=dpi,
    )

    if bar_groups and not top.empty:
        plot_top_variant_rpm_bars(
            df, sample_info,
            ids=top["ID"].astype(str).tolist(),
            groups=bar_groups,
            out_prefix=out_dir / f"RPM_bars_{tag}_top{len(top)}",
            title=f"Top {len(top)} {target} variants: RPM across groups",
            subtitle=(
                f"{subtitle}\nBars are the replicate median, whiskers the replicate "
                "range, points individual replicates; variants ordered by Er/Eh rank"
            ),
            dpi=dpi,
        )

    n_lib = int((~scores["is_control"]).sum())
    return {
        "arm": tag,
        "target": target,
        "backgrounds": "; ".join(backgrounds),
        "input_group": arm_input,
        "n_variants": n_lib,
        "n_in_box_median": int(
            ((~scores["is_control"])
             & (scores["Er_median"] > er_cut)
             & (scores["Eh_median"] > eh_cut)).sum()
        ),
        "n_in_box_mean": int(
            ((~scores["is_control"])
             & (scores["Er_mean"] > er_cut)
             & (scores["Eh_mean"] > eh_cut)).sum()
        ),
        "n_hits_median": int(len(hits)),
        "n_hits_mean": int(len(hits_mean)),
        "n_top_median": int(len(top)),
        "n_top_mean": int(len(top_mean)),
        "top_variant": top["ID"].iloc[0] if len(top) else np.nan,
        "top_Er_median": float(top["Er"].iloc[0]) if len(top) else np.nan,
        "top_Eh_median": float(top["Eh"].iloc[0]) if len(top) else np.nan,
    }


def run_enrichment_enhancement_analysis(
    df: pd.DataFrame,
    sample_info: pd.DataFrame,
    out_dir: Path,
    er_cut: float = 3.0,
    eh_cut: float = 3.0,
    min_rpm: float = 5.0,
    n_top: int = 10,
    pseudocount: float = 1.0,
    dpi: int = 300,
    write_excel: bool = False,
) -> pd.DataFrame:
    """Run enrichment/enhancement for all Sample_info arms and write a summary table."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    arms = build_enrichment_enhancement_arms(sample_info)
    write_table(arms, out_dir / "arms.tsv", write_excel=write_excel)
    summaries = []
    for _, arm in arms.iterrows():
        arm_input = str(arm["input_group"]) if "input_group" in arm.index else ""
        logger.info(
            "Enrichment/enhancement arm %s: target=%s backgrounds=%s input=%s",
            arm["tag"],
            arm["target"],
            arm["backgrounds"],
            arm_input,
        )
        summaries.append(
            run_enrichment_enhancement_arm(
                df, sample_info, arm, out_dir,
                input_group=arm_input,
                er_cut=er_cut, eh_cut=eh_cut, min_rpm=min_rpm,
                n_top=n_top, pseudocount=pseudocount, dpi=dpi,
                write_excel=write_excel,
            )
        )
    summary = pd.DataFrame(summaries)
    write_table(summary, out_dir / "Er_Eh_summary_all_arms.csv", write_excel=write_excel)
    return summary


def run_statistical_analysis(
    counts_file: Path,
    sample_info_file: Path,
    output_dir: Path,
    config_file: Optional[Path] = None,
    n_cpus: Optional[int] = None,
    verbose: bool = False,
    logger_instance: Optional[logging.Logger] = None,
    log_file: Optional[Path] = None,
) -> bool:
    """
    Run QC, DESeq2, and Lin-style enrichment/enhancement on CapScreen count data.

    QC files are written to ``output_dir/qc``. DESeq2 contrast tables are written
    to ``output_dir/deseq2``. Enrichment/enhancement results are written to
    ``output_dir/enrichment_enhancement``. The filtered count matrix, with RPM
    and CV columns, is written to ``output_dir/filtered_variant_metrics.tsv``.
    """
    if logger_instance is not None:
        global logger
        logger = logger_instance
        if log_file:
            _attach_stat_log_file(
                logger,
                Path(log_file),
                level=logging.DEBUG if verbose else logging.INFO,
            )
    else:
        setup_logging(verbose, log_file=log_file)
    
    config = load_config(config_file)
    output_dir.mkdir(parents=True, exist_ok=True)
    qc_dir = output_dir / "qc"
    de_dir = output_dir / "deseq2"
    ee_dir = output_dir / "enrichment_enhancement"
    qc_dir.mkdir(parents=True, exist_ok=True)
    de_dir.mkdir(parents=True, exist_ok=True)
    ee_dir.mkdir(parents=True, exist_ok=True)

    analysis_cfg = config.get("analysis", {}) or {}
    de_cfg = config.get("differential_expression", {}) or {}
    ee_cfg = config.get("enrichment_enhancement", {}) or {}
    output_cfg = config.get("output", {}) or {}
    plot_cfg = config.get("plots", {}) or {}
    n_cpus_val = n_cpus if n_cpus is not None else _cfg_int(analysis_cfg.get("n_cpus"), 8)
    min_count = de_cfg.get("min_count", analysis_cfg.get("low_expression_threshold", 10))
    min_count = 10 if min_count is None else float(min_count)
    drop_if_any = bool(de_cfg.get("drop_if_any_group_low", False))
    pseudocount = float(de_cfg.get("log2rpm_pseudocount", 1))
    min_replicates = _cfg_int(de_cfg.get("min_replicates"), 2)
    adjust_batch = bool(de_cfg.get("adjust_batch", True))
    write_excel = bool(output_cfg.get("write_excel", True))
    pca_top_n = _cfg_int((analysis_cfg.get("pca") or {}).get("top_n_features"), 500)
    ee_er_cut = float(ee_cfg.get("er_cut", 3))
    ee_eh_cut = float(ee_cfg.get("eh_cut", 3))
    ee_min_rpm = float(ee_cfg.get("min_rpm", 5))
    ee_n_top = _cfg_int(ee_cfg.get("n_top"), 10)
    ee_pseudo = float(ee_cfg.get("pseudocount", 1))
    plot_dpi = _cfg_int((plot_cfg.get("dpi") or {}).get("standard"), 300)
    pca_point_size = float((plot_cfg.get("pca") or {}).get("point_size", 120))
    pca_fontsize = float((plot_cfg.get("pca") or {}).get("fontsize", 10))
    pca_fig = (plot_cfg.get("figure_sizes") or {}).get("pca", [7.5, 6.5])
    if not isinstance(pca_fig, (list, tuple)) or len(pca_fig) != 2:
        pca_fig = [7.5, 6.5]
    pca_figsize = (float(pca_fig[0]), float(pca_fig[1]))

    analysis_params = {
        "analysis": {
            "n_cpus": n_cpus_val,
            "low_expression_threshold": min_count,
            "pca_top_n_features": pca_top_n,
        },
        "differential_expression": {
            "min_replicates": min_replicates,
            "min_count": min_count,
            "drop_if_any_group_low": drop_if_any,
            "log2rpm_pseudocount": pseudocount,
            "adjust_batch": adjust_batch,
        },
        "enrichment_enhancement": {
            "er_cut": ee_er_cut,
            "eh_cut": ee_eh_cut,
            "min_rpm": ee_min_rpm,
            "n_top": ee_n_top,
            "pseudocount": ee_pseudo,
        },
        "output": {
            "write_excel": write_excel,
        },
        "plots": {
            "dpi": plot_dpi,
            "pca_point_size": pca_point_size,
            "pca_fontsize": pca_fontsize,
            "pca_figsize": list(pca_figsize),
        },
    }
    with (output_dir / "analysis_params.json").open("w") as fh:
        json.dump(analysis_params, fh, indent=2)

    step = "initialization"
    try:
        step = "read count table"
        logger.info("STAT STEP START: %s", step)
        counts = load_merged_counts(counts_file)
        logger.info("STAT STEP OK: %s (%s)", step, counts.shape)

        step = "read sample info"
        logger.info("STAT STEP START: %s", step)
        meta = read_sample_info(sample_info_file)
        # Validate enrichment metadata early so DESeq2 and Er/Eh share one schema
        _validate_enrichment_metadata_consistency(meta)
        input_groups = list(
            dict.fromkeys(
                meta.loc[
                    meta["enrichment_role"].fillna("").astype(str).str.strip().str.lower().eq("input"),
                    "group",
                ]
                .astype(str)
                .tolist()
            )
        )
        logger.info(
            "STAT STEP OK: %s (%d samples, groups=%s, input=%s)",
            step,
            len(meta),
            sorted(meta["group"].astype(str).unique()),
            input_groups,
        )

        step = "collapse technical replicates"
        logger.info("STAT STEP START: %s", step)
        counts, meta = collapse_technical_replicates(counts, meta)
        _validate_enrichment_metadata_consistency(meta)
        input_groups = list(
            dict.fromkeys(
                meta.loc[
                    meta["enrichment_role"].fillna("").astype(str).str.strip().str.lower().eq("input"),
                    "group",
                ]
                .astype(str)
                .tolist()
            )
        )
        logger.info(
            "STAT STEP OK: %s (%d biological replicates, groups=%s)",
            step,
            len(meta),
            sorted(meta["group"].astype(str).unique()),
        )

        step = "variant metrics"
        logger.info("STAT STEP START: %s", step)
        filtered = update_variant_metrics(
            counts,
            meta,
            min_count=min_count,
            drop_if_any_group_low=drop_if_any,
            pseudocount=pseudocount,
        )
        filtered = add_cv_rpm_by_group(filtered, meta)
        filtered_path = output_dir / "filtered_variant_metrics.tsv"
        filtered.to_csv(filtered_path, sep="\t", index=False)
        logger.info("STAT STEP OK: %s (%s)", step, filtered_path)

        step = "log2RPM replicate correlation"
        logger.info("STAT STEP START: %s", step)
        try:
            corr = calc_pearson_biorep_by_group(filtered, meta)
            corr_path = qc_dir / "correlation_statistics.tsv"
            corr.to_csv(corr_path, sep="\t", index=False)
            plot_log2rpm_pearson_pdf(filtered, meta, qc_dir / "correlation_scatterPlot.pdf")
            logger.info("STAT STEP OK: %s (%s)", step, corr_path)
        except Exception as exc:
            logger.warning("STAT STEP FAILED: %s: %s (continuing)", step, _stat_failure_reason(step, exc))

        step = "CV histograms"
        logger.info("STAT STEP START: %s", step)
        try:
            plot_cv_histogram_pdf(filtered, qc_dir / "cv_rpm_histograms.pdf")
            logger.info("STAT STEP OK: %s", step)
        except Exception as exc:
            logger.warning("STAT STEP FAILED: %s: %s (continuing)", step, _stat_failure_reason(step, exc))

        step = "AAV2 contamination"
        logger.info("STAT STEP START: %s", step)
        try:
            plotted = plot_aav2_rpm_contamination(filtered, meta, qc_dir / "aav2_rpm_contamination.pdf")
            if plotted:
                logger.info("STAT STEP OK: %s", step)
            else:
                logger.warning("STAT STEP SKIPPED: %s", step)
        except Exception as exc:
            logger.warning("STAT STEP FAILED: %s: %s (continuing)", step, _stat_failure_reason(step, exc))

        step = "blind VST"
        logger.info("STAT STEP START: %s", step)
        vst, _meta_vst = fit_blind_vst(filtered, meta, n_cpus_val, adjust_batch=adjust_batch)
        logger.info("STAT STEP OK: %s (%s)", step, vst.shape)

        step = "VST Pearson heatmap"
        logger.info("STAT STEP START: %s", step)
        try:
            sample_cor = plot_vst_pearson_heatmap(vst, qc_dir / "vst_pearson_heatmap.pdf")
            sample_cor.to_csv(qc_dir / "vst_pearson_correlation.tsv", sep="\t")
            logger.info("STAT STEP OK: %s", step)
        except Exception as exc:
            logger.warning("STAT STEP FAILED: %s: %s (continuing)", step, _stat_failure_reason(step, exc))

        step = "VST PCA"
        logger.info("STAT STEP START: %s", step)
        try:
            plot_vst_pca(
                vst,
                meta,
                qc_dir,
                n_top=pca_top_n,
                exclude_groups=input_groups,
                point_size=pca_point_size,
                fontsize=pca_fontsize,
                figsize=pca_figsize,
            )
            logger.info("STAT STEP OK: %s", step)
        except Exception as exc:
            logger.warning("STAT STEP FAILED: %s: %s (continuing)", step, _stat_failure_reason(step, exc))

        step = "DESeq2 contrasts"
        logger.info("STAT STEP START: %s", step)
        try:
            comparisons = build_enrichment_comparisons(meta)
            write_table(comparisons, de_dir / "comparisons.tsv", write_excel=write_excel)
            n_ok = run_deseq2_comparisons(
                filtered,
                meta,
                comparisons,
                de_dir,
                n_cpus=n_cpus_val,
                min_replicates=min_replicates,
                adjust_batch=adjust_batch,
                write_excel=write_excel,
            )
            logger.info("STAT STEP OK: %s (%d/%d contrasts)", step, n_ok, len(comparisons))
        except Exception as exc:
            logger.warning("STAT STEP FAILED: %s: %s (continuing)", step, _stat_failure_reason(step, exc))

        step = "enrichment/enhancement"
        logger.info("STAT STEP START: %s (n_top=%d)", step, ee_n_top)
        try:
            summary = run_enrichment_enhancement_analysis(
                filtered,
                meta,
                ee_dir,
                er_cut=ee_er_cut,
                eh_cut=ee_eh_cut,
                min_rpm=ee_min_rpm,
                n_top=ee_n_top,
                pseudocount=ee_pseudo,
                dpi=plot_dpi,
                write_excel=write_excel,
            )
            logger.info(
                "STAT STEP OK: %s (%d arm(s); summary=%s)",
                step,
                len(summary),
                ee_dir / "Er_Eh_summary_all_arms.csv",
            )
        except Exception as exc:
            logger.warning("STAT STEP FAILED: %s: %s (continuing)", step, _stat_failure_reason(step, exc))

        logger.info("Statistical analysis completed successfully!")
        return True
    except Exception as exc:
        reason = _stat_failure_reason(step, exc)
        logger.error("STAT STEP FAILED: %s: %s", step, reason, exc_info=True)
        return False


def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(
        description="Statistical analysis of CapScreen count data",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--counts',
        type=Path,
        required=True,
        help='Path to merged.counts.csv file'
    )
    
    parser.add_argument(
        '--sample-info',
        type=Path,
        required=True,
        help='Path to Sample_info.csv file'
    )
    
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path.cwd(),
        help='Output directory for results (default: current directory)'
    )
    
    parser.add_argument(
        '--n-cpus',
        type=int,
        default=None,
        help='Number of CPUs to use for DESeq2 (default: from config.json or 8)'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    parser.add_argument(
        '--config',
        type=Path,
        default=None,
        help='Path to config.json file (default: looks for config.json in current directory)'
    )
    parser.add_argument(
        '--log-file',
        type=Path,
        default=None,
        help='Write stat step logs to this file (in addition to stdout)'
    )
    
    args = parser.parse_args()
    
    setup_logging(args.verbose, log_file=args.log_file)
    
    success = run_statistical_analysis(
        counts_file=args.counts,
        sample_info_file=args.sample_info,
        output_dir=args.output_dir,
        config_file=args.config,
        n_cpus=args.n_cpus,
        verbose=args.verbose,
        log_file=args.log_file,
    )
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()


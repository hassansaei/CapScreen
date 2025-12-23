#!/usr/bin/env python3
"""
Comprehensive statistical analysis script for CapScreen data.

This script performs a complete statistical analysis pipeline on CapScreen count data:

1. Data Processing:
   - Reads merged count tables (raw counts and RPM counts)
   - Processes sample metadata to identify groups, biological and technical replicates
   - Normalizes counts using DESeq2 (pydeseq2)

2. RPM-based Enrichment Analysis:
   - Computes mean RPM across technical replicates per biological replicate
   - Calculates log2 enrichment vs input groups
   - Aggregates enrichment statistics (mean, SD, CV) across experiments
   - Computes target specificity metrics

3. Correlation Analysis:
   - Calculates Spearman and Pearson correlations for all pairwise sample comparisons
   - Generates correlation scatter plots within each group

4. Principal Component Analysis (PCA):
   - PCA with all features
   - PCA with top variable features (configurable)
   - PCA excluding input samples

5. Differential Expression Analysis:
   - OLS regression-based differential expression (target vs background groups)
   - Calculates log2 fold changes, p-values, and FDR-adjusted p-values

6. Visualization:
   - Target enrichment scatter plots (with top enriched variants highlighted)
   - Top target-specific variants bar plots
   - RPM scatter density plots for replicate comparisons
   - Heatmaps for target-enriched peptides (from differential expression results)

All analysis parameters are configurable via config.json file.
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Tuple, Optional, Dict, Any

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from itertools import combinations
from scipy.stats import spearmanr, pearsonr
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import mpl_scatter_density  # noqa: F401

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
                "lambda": 0.01,
                "low_expression_threshold": 10,
                "pca": {
                    "n_components": 2,
                    "random_state": 0,
                    "top_n_features": 500
                }
            },
            "differential_expression": {
                "padj_threshold": 0.05,
                "log2fc_threshold": 1.0
            },
            "plots": {
                "scatter_highlight_top_n": 60,
                "bar_plot_top_n": 30,
                "dpi": {
                    "standard": 300,
                    "heatmap": 600
                },
                "figure_sizes": {
                    "scatter": [6, 6],
                    "bar": [6, 6],
                    "rpm_scatter": [6, 4],
                    "pca": [6, 5],
                    "heatmap": [10, 10]
                },
                "scatter": {
                    "all_points_size": 10,
                    "all_points_alpha": 0.4,
                    "highlighted_size": 40
                },
                "pca": {
                    "point_size": 120,
                    "fontsize": 10
                },
                "fontsize": 10
            },
            "colors": {
                "target": "#d62728",
                "background_bead": "#1f77b4",
                "background_bead_fc": "#ff7f0e",
                "gray": "#808080",
                "highlight": "#d62728"
            }
        }
    
    with open(config_file, 'r') as f:
        full_config = json.load(f)
    
    # Extract statistical_analysis section, or use full config if section doesn't exist
    config = full_config.get("statistical_analysis", full_config)
    
    logger.info(f"Loaded configuration from {config_file}")
    return config


def setup_logging(verbose: bool = False):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def read_count_table(count_file: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read merged count table and separate raw counts from RPM counts.
    
    Args:
        count_file: Path to merged.counts.tsv file
        
    Returns:
        Tuple of (raw_counts_df, rpm_counts_df)
    """
    logger.info(f"Reading count table from {count_file}")
    
    # Read the count table
    counts = pd.read_csv(count_file, sep=',')
    
    # Set index to ID_WLG
    if 'ID_WLG' in counts.columns:
        counts.index = counts['ID_WLG']
    else:
        logger.warning("ID_WLG column not found, using first column as index")
        counts.index = counts.iloc[:, 0]
    
    # Identify RPM columns
    rpm_columns = [col for col in counts.columns if col.endswith('_RPM')]
    
    # Create raw counts dataframe (drop RPM columns and metadata columns)
    raw_counts = counts.drop(columns=['ID_WLG', 'peptide'] + rpm_columns, errors='ignore')
    
    # Create RPM counts dataframe (keep only RPM columns)
    rpm_counts = counts[['ID_WLG', 'peptide'] + rpm_columns].copy() if rpm_columns else pd.DataFrame()
    
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
    
    meta = pd.read_csv(sample_info_file, sep=',', index_col=0)
    
    # Select relevant columns
    required_cols = ['group', 'tech_rep', 'bio_rep', 'batch_seq']
    missing_cols = [col for col in required_cols if col not in meta.columns]
    if missing_cols:
        logger.warning(f"Missing columns in sample info: {missing_cols}")
    
    meta = meta.loc[:, [col for col in required_cols if col in meta.columns]]
    
    logger.info(f"Sample info shape: {meta.shape}")
    logger.info(f"Groups found: {sorted(meta['group'].unique())}")
    
    return meta


def identify_groups(meta: pd.DataFrame) -> Tuple[list, list, list]:
    """
    Identify target, background, and input groups from metadata.
    
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
    # Background patterns: Bead, KO, Bead-Fc, Pool
    target_keywords = ['Target', 'target', 'WT', 'wt']
    background_keywords = ['Bead', 'bead', 'KO', 'ko', 'Pool', 'pool', 'Background', 'background']
    
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


def create_deseq_dataset(counts: pd.DataFrame, meta: pd.DataFrame, 
                        n_cpus: int = 8) -> DeseqDataSet:
    """
    Create and run DESeq2 dataset.
    
    Args:
        counts: Raw counts dataframe (samples as rows, variants as columns)
        meta: Metadata dataframe
        n_cpus: Number of CPUs to use
        
    Returns:
        DeseqDataSet object
    """
    logger.info("Creating DESeq2 dataset...")
    
    # Align counts and metadata
    common_samples = counts.index.intersection(meta.index)
    if len(common_samples) == 0:
        raise ValueError("No common samples found between counts and metadata")
    
    counts_aligned = counts.loc[common_samples]
    meta_aligned = meta.loc[common_samples]
    
    logger.info(f"Aligned {len(common_samples)} samples")
    
    # Create DESeqDataSet
    dds = DeseqDataSet(
        counts=counts_aligned,
        metadata=meta_aligned,
        design_factors=["group", "bio_rep"],  # batch-aware design
        refit_cooks=True,
        n_cpus=n_cpus
    )
    
    # Run DESeq2
    logger.info("Running DESeq2...")
    dds.deseq2()
    logger.info("DESeq2 completed")
    
    return dds


def extract_normalized_counts(dds: DeseqDataSet) -> pd.DataFrame:
    """
    Extract normalized counts from DESeqDataSet.
    
    Args:
        dds: DeseqDataSet object
        
    Returns:
        DataFrame with normalized counts (variants as rows, samples as columns)
    """
    logger.info("Extracting normalized counts...")
    
    norm_counts = dds.layers["normed_counts"]
    
    df_norm = pd.DataFrame(
        norm_counts,
        index=dds.obs_names,      # variant/peptide IDs (rows)
        columns=dds.var_names     # sample IDs (columns)
    )
    
    # Transpose to have variants as rows, samples as columns
    # (After transpose: rows=variants from dds.obs_names, cols=samples from dds.var_names)
    df_norm = df_norm.T
    
    # Remove AAV2 or other control variants if they exist in the index
    control_variants = df_norm.index[df_norm.index.str.contains('AAV', case=False, na=False)]
    if len(control_variants) > 0:
        logger.info(f"Removing control variants: {list(control_variants)}")
        df_norm = df_norm.drop(control_variants)
    
    logger.info(f"Normalized counts shape: {df_norm.shape}")
    
    return df_norm


def calculate_correlation_statistics(dds: DeseqDataSet, df_norm: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate correlation statistics (Spearman and Pearson) for each group.
    
    Args:
        dds: DeseqDataSet object
        df_norm: Normalized counts dataframe (variants as rows, samples as columns)
        
    Returns:
        DataFrame with correlation statistics for each group
    """
    logger.info("Calculating correlation statistics...")
    
    corr_results = {}
    groups = dds.obs['group'].unique()
    
    for g in groups:
        samples = dds.obs.index[dds.obs['group'] == g].tolist()
        
        if len(samples) >= 2:  # at least 2 reps needed
            pair_s = []
            pair_p = []
            
            for s1, s2 in combinations(samples, 2):
                # Ensure both samples exist in df_norm (samples are columns after transpose)
                if s1 in df_norm.columns and s2 in df_norm.columns:
                    s_val, _ = spearmanr(df_norm[s1], df_norm[s2])
                    p_val, _ = pearsonr(df_norm[s1], df_norm[s2])
                    pair_s.append(s_val)
                    pair_p.append(p_val)
            
            if pair_s:  # Only add if we have correlations
                corr_results[g] = {
                    "spearman_mean": np.mean(pair_s),
                    "spearman_min": np.min(pair_s),
                    "spearman_max": np.max(pair_s),
                    "spearman_sd": np.std(pair_s),
                    "pearson_mean": np.mean(pair_p),
                    "pearson_min": np.min(pair_p),
                    "pearson_max": np.max(pair_p),
                    "pearson_sd": np.std(pair_p)
                }
    
    df_corr = pd.DataFrame.from_dict(corr_results, orient='index')
    
    logger.info(f"Calculated correlations for {len(corr_results)} groups")
    
    return df_corr


def create_correlation_scatter_plots(dds: DeseqDataSet, df_norm: pd.DataFrame, 
                                     output_file: Path):
    """
    Create scatter plots for all pairwise sample comparisons within each group.
    
    Args:
        dds: DeseqDataSet object
        df_norm: Normalized counts dataframe (variants as rows, samples as columns)
        output_file: Path to save the PDF file
    """
    logger.info(f"Creating correlation scatter plots...")
    
    groups = dds.obs['group'].unique()
    
    with PdfPages(output_file) as pdf:
        for g in groups:
            samples = dds.obs.index[dds.obs['group'] == g].tolist()
            
            if len(samples) >= 2:
                for s1, s2 in combinations(samples, 2):
                    # Ensure both samples exist in df_norm (samples are columns after transpose)
                    if s1 in df_norm.columns and s2 in df_norm.columns:
                        sp, _ = spearmanr(df_norm[s1], df_norm[s2])
                        pr, _ = pearsonr(df_norm[s1], df_norm[s2])
                        
                        fig, ax = plt.subplots()
                        ax.scatter(df_norm[s1], df_norm[s2], color="darkblue", s=15)
                        ax.set_title(f"{g} | Spearman r={sp:.2f} | Pearson r={pr:.2f}")
                        ax.set_xlabel(s1)
                        ax.set_ylabel(s2)
                        fig.tight_layout()
                        pdf.savefig(fig)
                        plt.close(fig)
    
    logger.info(f"Saved correlation scatter plots to {output_file}")


def perform_pca(df_norm: pd.DataFrame, meta: pd.DataFrame, 
                top_n_features: Optional[int] = None,
                exclude_input: bool = False,
                input_groups: Optional[list] = None,
                config: Optional[Dict[str, Any]] = None) -> Tuple[PCA, pd.DataFrame]:
    """
    Perform PCA analysis on normalized counts.
    
    Args:
        df_norm: Normalized counts dataframe (variants as rows, samples as columns)
        meta: Metadata dataframe
        top_n_features: Number of top variable features to use (None = all features)
        exclude_input: Whether to exclude input samples
        input_groups: List of input group names (used if exclude_input=True)
        
    Returns:
        Tuple of (PCA object, DataFrame with PC coordinates and metadata)
    """
    # Align metadata (samples are columns in df_norm, so we need to align with columns)
    common_samples = df_norm.columns.intersection(meta.index)
    df_norm_aligned = df_norm[common_samples]  # Select columns (samples)
    meta_aligned = meta.loc[common_samples]
    
    # Filter out input samples if requested (before filtering variants)
    if exclude_input and input_groups:
        samples_to_keep = meta_aligned[~meta_aligned['group'].isin(input_groups)].index
        df_norm_aligned = df_norm_aligned[samples_to_keep]
        meta_aligned = meta_aligned.loc[samples_to_keep]
        logger.info(f"Excluded input samples. Remaining samples: {len(samples_to_keep)}")
    
    # Get config values
    if config is None:
        config = {}
    low_expr_threshold = config.get("analysis", {}).get("low_expression_threshold", 10)
    pca_config = config.get("analysis", {}).get("pca", {})
    n_components = pca_config.get("n_components", 2)
    random_state = pca_config.get("random_state", 0)
    
    # Filter low-enriched peptides (variants as rows, sum across samples)
    df_norm_aligned = df_norm_aligned.loc[df_norm_aligned.sum(axis=1) > low_expr_threshold]
    
    # Log transform
    log_counts = np.log2(df_norm_aligned + 1)
    
    # Select top variable features if specified (variants are rows, so axis=1)
    if top_n_features is not None:
        feature_var = log_counts.var(axis=1, ddof=1)  # Variance across samples for each variant
        top_features = feature_var.nlargest(top_n_features).index
        log_counts = log_counts.loc[top_features]
        logger.info(f"Using top {top_n_features} variable features")
    
    # PCA: transpose so samples are rows, variants are columns
    X = log_counts.T  # samples x features
    X_scaled = StandardScaler(with_mean=True, with_std=True).fit_transform(X)
    
    # Perform PCA
    pca = PCA(n_components=n_components, random_state=random_state)
    pcs = pca.fit_transform(X_scaled)
    
    # Create results dataframe
    pca_df = pd.DataFrame(
        pcs,
        index=X.index,  # sample names
        columns=["PC1", "PC2"]
    )
    
    # Merge with metadata
    pca_df = pca_df.join(meta_aligned)
    
    logger.info(f"PCA variance explained: PC1={pca.explained_variance_ratio_[0]:.3f}, "
                f"PC2={pca.explained_variance_ratio_[1]:.3f}")
    
    return pca, pca_df


def perform_differential_expression(df_norm: pd.DataFrame, meta: pd.DataFrame,
                                    target_group: str, background_group: str) -> pd.DataFrame:
    """
    Perform differential expression analysis using OLS regression.
    
    Args:
        df_norm: Normalized counts dataframe (variants as rows, samples as columns)
        meta: Metadata dataframe
        target_group: Name of target group
        background_group: Name of background group
        
    Returns:
        DataFrame with differential expression results (log2FC, pval, padj)
    """
    logger.info(f"Performing differential expression: {target_group} vs {background_group}")
    
    # 1) Subset to target vs background groups
    groups_de = [target_group, background_group]
    meta_de = meta[meta["group"].isin(groups_de)].copy()
    samples_de = meta_de.index
    
    if len(samples_de) == 0:
        raise ValueError(f"No samples found for groups {groups_de}")
    
    logger.info(f"Found {len(samples_de)} samples for comparison")
    
    # 2) Build response matrix (ensure numeric + finite)
    # df_norm has variants as rows, samples as columns
    Y = df_norm[samples_de].apply(pd.to_numeric, errors="coerce")
    Y = np.log2(Y + 1)
    
    # Drop peptides with any NaN/Inf
    Y = Y.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    
    logger.info(f"Using {len(Y)} peptides for analysis")
    
    # 3) Design matrix (ensure numeric + aligned)
    X = pd.DataFrame(index=samples_de)
    X[target_group] = (meta_de.loc[samples_de, "group"] == target_group).astype(int).astype(float)
    
    # Add biological replicate dummies (drop first to avoid multicollinearity)
    bio = pd.get_dummies(meta_de.loc[samples_de, "bio_rep"].astype(str), 
                         prefix="bio", drop_first=True)
    bio = bio.astype(float)
    
    X = pd.concat([X, bio], axis=1)
    X = sm.add_constant(X).astype(float)
    
    # Final sanity checks
    assert np.isfinite(X.to_numpy()).all(), "Design matrix contains non-finite values"
    assert np.isfinite(Y.to_numpy()).all(), "Response matrix contains non-finite values"
    
    # 4) Fit per peptide
    results = []
    X_np = X.to_numpy()
    target_col = list(X.columns).index(target_group)
    
    logger.info("Fitting OLS models for each peptide...")
    for pep in Y.index:
        y = Y.loc[pep].to_numpy(dtype=float)
        try:
            fit = sm.OLS(y, X_np, missing="drop").fit()
            coef = fit.params[target_col]
            pval = fit.pvalues[target_col]
            results.append((pep, coef, pval))
        except Exception as e:
            logger.warning(f"Failed to fit model for {pep}: {e}")
            continue
    
    # Create results dataframe
    res = pd.DataFrame(results, columns=["peptide", f"log2FC_{target_group}_vs_{background_group}", "pval"])
    res = res.set_index("peptide")
    
    # Calculate adjusted p-values
    res["padj"] = multipletests(res["pval"].values, method="fdr_bh")[1]
    
    # Sort by adjusted p-value
    res_sorted = res.sort_values("padj")
    
    logger.info(f"Completed differential expression analysis. "
                f"Found {len(res_sorted[res_sorted['padj'] < 0.05])} significant peptides (padj < 0.05)")
    
    return res_sorted


def compute_rpm_enrichment_analysis(rpm_counts: pd.DataFrame, meta: pd.DataFrame,
                                    input_groups: list, LAMBDA: float = 0.01) -> pd.DataFrame:
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
        
    Returns:
        Annotated RPM dataframe with enrichment metrics
    """
    logger.info("Computing RPM-based enrichment analysis...")
    
    # Prepare RPM counts: set ID_WLG as index if it's a column
    df_rpm = rpm_counts.copy()
    if 'ID_WLG' in df_rpm.columns:
        df_rpm = df_rpm.set_index('ID_WLG')
        df_rpm.index.name = 'peptide'  # Rename index for consistency
    elif df_rpm.index.name != 'peptide':
        df_rpm.index.name = 'peptide'
    
    if 'peptide' in df_rpm.columns:
        df_rpm = df_rpm.drop(columns='peptide', errors='ignore')
    
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
            mu.append(
                pd.DataFrame({
                    "peptide": mu_r.index,
                    "mu": mu_r.values,
                    "group": group,
                    "bio_rep": rep
                })
            )
    
    mu_df = pd.concat(mu, ignore_index=True)
    logger.info(f"Computed means for {len(mu_df)} peptide-group-bio_rep combinations")
    
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
        
        # Average across input groups if multiple exist
        mu_input = mu_input_df.groupby("peptide")["mu"].mean()
        
        for group in non_input_groups:
            mu_g = (
                mu_r[mu_r["group"] == group]
                .set_index("peptide")["mu"]
            )
            
            common = mu_g.index.intersection(mu_input.index)
            
            if len(common) > 0:
                e = np.log2((mu_g.loc[common] + LAMBDA) / (mu_input.loc[common] + LAMBDA))
                
                enrich.append(
                    pd.DataFrame({
                        "peptide": common,
                        "log2_enrichment": e.values,
                        "group": group,
                        "bio_rep": rep
                    })
                )
    
    if len(enrich) == 0:
        logger.warning("No enrichment data computed. Check input groups.")
        return df_rpm
    
    enrich_df = pd.concat(enrich, ignore_index=True)
    logger.info(f"Computed enrichments for {len(enrich_df)} peptide-group-bio_rep combinations")
    
    # Step 3: Aggregate across experiments
    logger.info("Step 3: Aggregating across experiments...")
    enrich_summary = (
        enrich_df
        .groupby(["peptide", "group"])
        .agg(
            mean_log2_enrichment=("log2_enrichment", "mean"),
            sd_log2_enrichment=("log2_enrichment", "std")
        )
        .reset_index()
    )
    
    # Step 4: Create wide format tables and compute statistics
    logger.info("Step 4: Creating summary statistics...")
    df_rpm_mat = df_rpm.copy()
    
    # Pivot mean and sd
    mean_wide = enrich_summary.pivot(index="peptide", columns="group", values="mean_log2_enrichment")
    sd_wide = enrich_summary.pivot(index="peptide", columns="group", values="sd_log2_enrichment")
    
    # Rename columns
    # Calculate CV (coefficient of variation) with protection against division by zero
    cv_wide = sd_wide / mean_wide.abs().replace(0, np.nan)
    cv_wide = cv_wide.fillna(0)  # Set CV to 0 where mean is 0
    cv_wide = cv_wide.add_prefix("CV_")
    mean_wide = mean_wide.add_prefix("mean_log2Enrich_")
    sd_wide = sd_wide.add_prefix("sd_log2Enrich_")
    
    # Combine
    peptide_annot = pd.concat([mean_wide, sd_wide, cv_wide], axis=1)
    
    # Compute target specificity if we have Target and Bead-Fc groups
    target_groups = [g for g in mean_wide.columns if 'Target' in g or 'target' in g]
    background_groups = [g for g in mean_wide.columns if any(x in g for x in ['Bead', 'bead', 'KO', 'ko'])]
    
    if target_groups and background_groups:
        target_col = f"mean_log2Enrich_{target_groups[0]}"
        background_col = f"mean_log2Enrich_{background_groups[0]}"
        if target_col in peptide_annot.columns and background_col in peptide_annot.columns:
            # Create specificity column name (handle special characters like "-")
            target_name = target_groups[0].replace("-", "").replace("_", "")
            bg_name = background_groups[0].replace("-", "").replace("_", "")
            specificity_col_name = f"{target_name}_minus_{bg_name}"
            peptide_annot[specificity_col_name] = (
                peptide_annot[target_col] - peptide_annot[background_col]
            )
            logger.info(f"Created specificity column: {specificity_col_name}")
    
    # Reset index to merge (peptide_annot index is already "peptide")
    peptide_annot = peptide_annot.reset_index()
    
    # Merge into df_rpm (df_rpm_mat index should be peptide IDs)
    df_rpm_annot = (
        df_rpm_mat
        .reset_index()  # Reset to make index a column for merging
        .merge(peptide_annot, on="peptide", how="left")
        .set_index("peptide")
    )
    
    # Step 5: Add per-replicate enrichments in wide format
    logger.info("Step 5: Adding per-replicate enrichments...")
    enrich_wide = (
        enrich_df
        .pivot_table(
            index="peptide",
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


def plot_target_enrichment_scatter(df_rpm_annot: pd.DataFrame, output_file: Path,
                                   target_group: str = "Target", 
                                   background_group: str = "Bead-Fc",
                                   top_n: int = 60,
                                   config: Optional[Dict[str, Any]] = None):
    """
    Create scatter plot of target vs background enrichment, highlighting top enriched variants.
    
    Args:
        df_rpm_annot: Annotated RPM dataframe with enrichment metrics
        output_file: Path to save the plot
        target_group: Name of target group (default: "Target")
        background_group: Name of background group (default: "Bead-Fc")
        top_n: Number of top peptides to highlight (default: 60)
    """
    logger.info(f"Creating target enrichment scatter plot: {target_group} vs {background_group}")
    
    # Column names for enrichment
    target_col = f"mean_log2Enrich_{target_group}"
    background_col = f"mean_log2Enrich_{background_group}"
    
    # Check if columns exist
    if target_col not in df_rpm_annot.columns or background_col not in df_rpm_annot.columns:
        logger.warning(f"Required columns not found: {target_col} or {background_col}")
        logger.warning(f"Available columns: {list(df_rpm_annot.columns)}")
        return
    
    # Filter out NaN values and compute specificity on the fly
    df_plot = df_rpm_annot[[target_col, background_col]].dropna()
    
    if len(df_plot) == 0:
        logger.warning("No data available for scatter plot")
        return
    
    # Compute specificity: target - background
    df_plot = df_plot.copy()
    df_plot['specificity'] = df_plot[target_col] - df_plot[background_col]
    
    # Get config values
    if config is None:
        config = {}
    plot_config = config.get("plots", {})
    figsize = plot_config.get("figure_sizes", {}).get("scatter", [6, 6])
    dpi = plot_config.get("dpi", {}).get("standard", 300)
    scatter_config = plot_config.get("scatter", {})
    all_size = scatter_config.get("all_points_size", 10)
    all_alpha = scatter_config.get("all_points_alpha", 0.4)
    highlight_size = scatter_config.get("highlighted_size", 40)
    fontsize = plot_config.get("fontsize", 10)
    colors_config = config.get("colors", {})
    gray_color = colors_config.get("gray", "#808080")
    highlight_color = colors_config.get("highlight", "#d62728")
    
    # Create figure
    plt.figure(figsize=tuple(figsize), dpi=dpi)
    
    # Plot all points in gray
    plt.scatter(
        df_plot[background_col],
        df_plot[target_col],
        s=all_size, alpha=all_alpha, color=gray_color
    )
    
    # Highlight top peptides by sorting on specificity (target - background)
    top = df_plot.sort_values('specificity', ascending=False).head(top_n)
    
    if len(top) > 0:
        plt.scatter(
            top[background_col],
            top[target_col],
            s=highlight_size, color=highlight_color, zorder=5  # zorder to ensure points are on top
        )
        logger.info(f"Highlighted {len(top)} top enriched variants (sorted by {target_col} - {background_col})")
    else:
        logger.warning(f"No peptides found for highlighting")
    
    # Add diagonal line
    plt.axline((0, 0), slope=1, linestyle="--", color="black")
    
    # Set labels
    plt.xlabel(f"{background_group} log2 enrichment vs Input", fontsize=fontsize)
    plt.ylabel(f"{target_group} log2 enrichment vs Input", fontsize=fontsize)
    plt.title("Target-specific peptide enrichment", fontsize=fontsize)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved target enrichment scatter plot to {output_file}")


def plot_top_target_specific_variants(df_rpm_annot: pd.DataFrame, output_file: Path,
                                      target_group: str = "Target",
                                      background_group: str = "Bead-Fc",
                                      top_n: int = 30,
                                      config: Optional[Dict[str, Any]] = None):
    """
    Create horizontal bar plot showing top target-specific variants.
    
    Args:
        df_rpm_annot: Annotated RPM dataframe with enrichment metrics
        output_file: Path to save the plot
        target_group: Name of target group (default: "Target")
        background_group: Name of background group (default: "Bead-Fc")
        top_n: Number of top variants to show (default: 30)
    """
    logger.info(f"Creating top {top_n} target-specific variants bar plot: {target_group} vs {background_group}")
    
    # Column names for enrichment
    target_col = f"mean_log2Enrich_{target_group}"
    background_col = f"mean_log2Enrich_{background_group}"
    
    # Check if columns exist
    if target_col not in df_rpm_annot.columns or background_col not in df_rpm_annot.columns:
        logger.warning(f"Required columns not found: {target_col} or {background_col}")
        logger.warning(f"Available columns: {list(df_rpm_annot.columns)}")
        return
    
    # Compute specificity column name (matching the pattern from compute_rpm_enrichment_analysis)
    bg_name_clean = background_group.replace("-", "").replace("_", "")
    target_name_clean = target_group.replace("-", "").replace("_", "")
    specificity_col_name = f"{target_name_clean}_minus_{bg_name_clean}"
    
    # Check if specificity column exists, if not compute it
    if specificity_col_name not in df_rpm_annot.columns:
        logger.info(f"Computing specificity column: {specificity_col_name}")
        df_plot = df_rpm_annot[[target_col, background_col]].dropna().copy()
        df_plot[specificity_col_name] = df_plot[target_col] - df_plot[background_col]
    else:
        logger.info(f"Using existing specificity column: {specificity_col_name}")
        df_plot = df_rpm_annot[[target_col, background_col, specificity_col_name]].dropna().copy()
    
    if len(df_plot) == 0:
        logger.warning("No data available for bar plot")
        return
    
    # Sort by specificity and get top N
    top_variants = df_plot.sort_values(specificity_col_name, ascending=False).head(top_n)
    
    if len(top_variants) == 0:
        logger.warning("No variants found for bar plot")
        return
    
    # Get config values
    if config is None:
        config = {}
    plot_config = config.get("plots", {})
    figsize = plot_config.get("figure_sizes", {}).get("bar", [6, 6])
    dpi = plot_config.get("dpi", {}).get("standard", 300)
    fontsize = plot_config.get("fontsize", 10)
    colors_config = config.get("colors", {})
    highlight_color = colors_config.get("highlight", "#d62728")
    
    # Create figure
    plt.figure(figsize=tuple(figsize), dpi=dpi)
    
    # Create horizontal bar plot
    plt.barh(
        top_variants.index,
        top_variants[specificity_col_name],
        color=highlight_color
    )
    
    # Set labels
    plt.xlabel(f"Target − {background_group} log2 enrichment", fontsize=fontsize)
    plt.ylabel("Peptide", fontsize=fontsize)
    plt.title(f"Top {len(top_variants)} Target-specific peptides", fontsize=fontsize)
    
    # Invert y-axis so highest values are at the top
    plt.gca().invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved top {len(top_variants)} target-specific variants bar plot to {output_file}")


def create_white_viridis_colormap() -> LinearSegmentedColormap:
    """
    Create a custom colormap with white background for scatter density plots.
    
    Returns:
        LinearSegmentedColormap object
    """
    white_viridis = LinearSegmentedColormap.from_list(
        "white_viridis",
        [
            (0, "#ffffff"),
            (1e-20, "#440053"),
            (0.2, "#404388"),
            (0.4, "#2a788e"),
            (0.6, "#21a784"),
            (0.8, "#78d151"),
            (1, "#fde624"),
        ],
        N=256,
    )
    return white_viridis


def create_rpm_scatter_density_plots(rpm_counts: pd.DataFrame, meta: pd.DataFrame,
                                     output_dir: Path,
                                     config: Optional[Dict[str, Any]] = None):
    """
    Create scatter density plots for all replicate comparisons within each group using RPM counts.
    
    Args:
        rpm_counts: RPM counts dataframe (variants as rows, samples as columns)
        meta: Metadata dataframe
        output_dir: Directory to save PNG plots
    """
    logger.info("Creating scatter density plots for RPM counts...")
    
    # Reset matplotlib style to default with white background
    plt.style.use("default")
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })
    
    # Create custom colormap
    white_viridis = create_white_viridis_colormap()
    
    # Align RPM counts with metadata
    # rpm_counts has variants as rows, samples as columns
    # Remove ID_WLG and peptide columns if they exist
    if 'ID_WLG' in rpm_counts.columns:
        rpm_counts = rpm_counts.set_index('ID_WLG')
    if 'peptide' in rpm_counts.columns:
        rpm_counts = rpm_counts.drop(columns='peptide', errors='ignore')
    
    # Remove AAV2 or control variants if they exist
    control_variants = rpm_counts.index[rpm_counts.index.str.contains('AAV', case=False, na=False)]
    if len(control_variants) > 0:
        logger.info(f"Removing control variants from RPM counts: {list(control_variants)}")
        rpm_counts = rpm_counts.drop(control_variants)
    
    # Get common samples
    common_samples = rpm_counts.columns.intersection(meta.index)
    rpm_counts_aligned = rpm_counts[common_samples]
    meta_aligned = meta.loc[common_samples]
    
    logger.info(f"Aligned {len(common_samples)} samples for RPM scatter density plots")
    
    # Get all groups
    groups = meta_aligned['group'].unique()
    
    plot_count = 0
    for group in groups:
        group_samples = meta_aligned[meta_aligned['group'] == group].index.tolist()
        
        if len(group_samples) >= 2:
            # Create plots for all pairwise comparisons within the group
            for s1, s2 in combinations(group_samples, 2):
                if s1 in rpm_counts_aligned.columns and s2 in rpm_counts_aligned.columns:
                    # Get RPM values and apply log2 transformation
                    x_rpm = rpm_counts_aligned[s1].to_numpy()
                    y_rpm = rpm_counts_aligned[s2].to_numpy()
                    
                    # Apply log2 transformation (add 1 to avoid log(0))
                    x = np.log2(x_rpm + 1)
                    y = np.log2(y_rpm + 1)
                    
                    # Remove NaN and Inf values
                    mask = np.isfinite(x) & np.isfinite(y)
                    x = x[mask]
                    y = y[mask]
                    
                    if len(x) > 0:
                        # Get config values
                        if config is None:
                            config = {}
                        plot_config = config.get("plots", {})
                        figsize = plot_config.get("figure_sizes", {}).get("rpm_scatter", [6, 4])
                        dpi = plot_config.get("dpi", {}).get("standard", 300)
                        fontsize = plot_config.get("fontsize", 10)
                        
                        fig = plt.figure(figsize=tuple(figsize), dpi=dpi)
                        ax = fig.add_subplot(1, 1, 1, projection="scatter_density")
                        
                        # Force white axes background
                        ax.patch.set_facecolor("white")
                        
                        # Create scatter density plot
                        density = ax.scatter_density(x, y, cmap=white_viridis)
                        
                        # Add colorbar
                        cbar = fig.colorbar(density, ax=ax)
                        cbar.set_label("Density")
                        
                        # Set labels and title
                        ax.set_xlabel(s1, fontsize=fontsize)
                        ax.set_ylabel(s2, fontsize=fontsize)
                        ax.set_title(f"{group}", fontsize=fontsize)
                        
                        plt.tight_layout()
                        
                        # Save plot
                        safe_group_name = group.replace(" ", "_").replace("/", "_")
                        plot_filename = f"rpm_scatter_{safe_group_name}_{s1}_vs_{s2}.png"
                        plot_path = output_dir / plot_filename
                        plt.savefig(plot_path, dpi=dpi, bbox_inches='tight')
                        plt.close()
                        
                        plot_count += 1
    
    logger.info(f"Created {plot_count} scatter density plots for RPM counts")


def create_group_color_palette(all_groups: list) -> dict:
    """
    Create a fixed color palette for groups to ensure consistency across plots.
    
    Args:
        all_groups: List of all unique group names
        
    Returns:
        Dictionary mapping group names to colors
    """
    # Use a consistent color palette (tab10 for up to 10 groups, then extend)
    base_colors = sns.color_palette("tab10", n_colors=10)
    
    # If more than 10 groups, extend with additional colors
    if len(all_groups) > 10:
        extended_colors = sns.color_palette("Set3", n_colors=len(all_groups) - 10)
        colors = base_colors + extended_colors
    else:
        colors = base_colors[:len(all_groups)]
    
    # Create mapping: sort groups for consistency
    sorted_groups = sorted(all_groups)
    palette = {group: colors[i] for i, group in enumerate(sorted_groups)}
    
    return palette


def plot_pca(pca: PCA, pca_df: pd.DataFrame, output_file: Path, 
             title: str = "PCA Plot", group_palette: Optional[dict] = None,
             config: Optional[Dict[str, Any]] = None):
    """
    Create and save PCA plot.
    
    Args:
        pca: PCA object
        pca_df: DataFrame with PC coordinates and metadata
        output_file: Path to save the plot
        title: Plot title
        group_palette: Dictionary mapping group names to colors (for consistency)
    """
    logger.info(f"Creating PCA plot: {title}")
    
    # Get config values
    if config is None:
        config = {}
    plot_config = config.get("plots", {})
    figsize = plot_config.get("figure_sizes", {}).get("pca", [6, 5])
    dpi = plot_config.get("dpi", {}).get("standard", 300)
    pca_plot_config = plot_config.get("pca", {})
    point_size = pca_plot_config.get("point_size", 120)
    fontsize = pca_plot_config.get("fontsize", 10)
    
    # PCA plot
    plt.figure(figsize=tuple(figsize), dpi=dpi)
    
    # Use fixed palette if provided
    plot_kwargs = {
        "data": pca_df,
        "x": "PC1",
        "y": "PC2",
        "hue": "group",
        "style": "bio_rep",
        "s": point_size
    }
    
    if group_palette is not None:
        plot_kwargs["palette"] = group_palette
    
    sns.scatterplot(**plot_kwargs)
    
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)", fontsize=fontsize)
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)", fontsize=fontsize)
    plt.legend(title="", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved PCA plot to {output_file}")


def run_statistical_analysis(
    counts_file: Path,
    sample_info_file: Path,
    output_dir: Path,
    config_file: Optional[Path] = None,
    n_cpus: Optional[int] = None,
    verbose: bool = False,
    logger_instance: Optional[logging.Logger] = None
) -> bool:
    """
    Run statistical analysis on CapScreen count data.
    
    This function can be called programmatically from other modules.
    
    Args:
        counts_file: Path to merged.counts.tsv file
        sample_info_file: Path to Sample_info.csv file
        output_dir: Output directory for results
        config_file: Path to config.json file (optional)
        n_cpus: Number of CPUs to use for DESeq2 (optional)
        verbose: Enable verbose logging (default: False)
        logger_instance: Optional logger instance to use (if None, uses module logger)
        
    Returns:
        True if analysis completed successfully, False otherwise
    """
    # Use provided logger or module logger
    if logger_instance is not None:
        global logger
        logger = logger_instance
    
    # Set up logging if not already configured
    if not logger.handlers:
        setup_logging(verbose)
    
    # Load configuration
    config = load_config(config_file)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Read data
        raw_counts, rpm_counts = read_count_table(counts_file)
        meta = read_sample_info(sample_info_file)
        
        # Identify groups
        target_groups, background_groups, input_groups = identify_groups(meta)
        
        # Determine target and background groups for differential expression
        # Use first target and first background group, or infer from group names
        if len(target_groups) > 0 and len(background_groups) > 0:
            de_target = target_groups[0]
            de_background = background_groups[0]
        elif len(target_groups) >= 2:
            # If we have multiple non-input groups, use first two
            de_target = target_groups[0]
            de_background = target_groups[1] if len(target_groups) > 1 else None
        else:
            # Try to infer from group names (common patterns)
            all_non_input = target_groups
            # Look for common patterns: Target, WT, etc. for target; Bead, KO, etc. for background
            target_candidates = [g for g in all_non_input if any(x in g for x in ['Target', 'WT', 'target'])]
            background_candidates = [g for g in all_non_input if any(x in g for x in ['Bead', 'KO', 'background', 'Background'])]
            
            if target_candidates and background_candidates:
                de_target = target_candidates[0]
                de_background = background_candidates[0]
            elif len(all_non_input) >= 2:
                de_target = all_non_input[0]
                de_background = all_non_input[1]
            else:
                logger.warning("Could not determine target and background groups for differential expression")
                de_target = None
                de_background = None
        
        # Create DESeq2 dataset
        n_cpus_val = n_cpus if n_cpus is not None else config.get("analysis", {}).get("n_cpus", 8)
        dds = create_deseq_dataset(raw_counts, meta, n_cpus=n_cpus_val)
        
        # Extract normalized counts (AAV2 control already removed in extract_normalized_counts)
        df_norm = extract_normalized_counts(dds)
        
        # Save normalized counts
        norm_counts_file = output_dir / "normalized_counts.tsv"
        df_norm.to_csv(norm_counts_file, sep='\t')
        logger.info(f"Saved normalized counts to {norm_counts_file}")
        
        # Create RPM scatter density plots for replicate comparisons
        if not rpm_counts.empty:
            logger.info("Creating scatter density plots for RPM counts...")
            try:
                create_rpm_scatter_density_plots(rpm_counts, meta, output_dir, config=config)
            except Exception as e:
                logger.warning(f"Failed to create RPM scatter density plots: {e}")
        else:
            logger.warning("RPM counts not available, skipping scatter density plots")
        
        # Perform RPM-based enrichment analysis
        df_rpm_annot = None
        if not rpm_counts.empty and input_groups:
            logger.info("Performing RPM-based enrichment analysis...")
            try:
                lambda_val = config.get("analysis", {}).get("lambda", 0.01)
                df_rpm_annot = compute_rpm_enrichment_analysis(
                    rpm_counts, meta, input_groups, LAMBDA=lambda_val
                )
                
                # Save annotated RPM table with informative name
                rpm_enrichment_file = output_dir / "RPM_enrichment_analysis.tsv"
                df_rpm_annot.to_csv(rpm_enrichment_file, sep='\t')
                logger.info(f"Saved RPM enrichment analysis to {rpm_enrichment_file}")
                
                # Create target enrichment scatter plot
                if df_rpm_annot is not None:
                    # Try to find target and background groups from the enrichment columns
                    enrich_cols = [col for col in df_rpm_annot.columns if col.startswith("mean_log2Enrich_")]
                    if enrich_cols:
                        # Extract group names from column names
                        groups_in_enrich = [col.replace("mean_log2Enrich_", "") for col in enrich_cols]
                        
                        # Find target and background groups
                        plot_target = None
                        plot_background = None
                        
                        for g in groups_in_enrich:
                            if any(x in g for x in ['Target', 'target', 'WT', 'wt']):
                                plot_target = g
                            elif any(x in g for x in ['Bead', 'bead', 'KO', 'ko', 'Background', 'background']):
                                plot_background = g
                        
                        if plot_target and plot_background:
                            scatter_file = output_dir / "target_enrichment_scatter.png"
                            scatter_top_n = config.get("plots", {}).get("scatter_highlight_top_n", 60)
                            plot_target_enrichment_scatter(
                                df_rpm_annot, scatter_file,
                                target_group=plot_target,
                                background_group=plot_background,
                                top_n=scatter_top_n,
                                config=config
                            )
                            
                            # Create top variants bar plot
                            bar_plot_file = output_dir / "top_target_specific_variants.png"
                            bar_top_n = config.get("plots", {}).get("bar_plot_top_n", 30)
                            plot_top_target_specific_variants(
                                df_rpm_annot, bar_plot_file,
                                target_group=plot_target,
                                background_group=plot_background,
                                top_n=bar_top_n,
                                config=config
                            )
                        else:
                            logger.warning("Could not determine target/background groups for scatter plot")
                
            except Exception as e:
                logger.warning(f"Failed to perform RPM enrichment analysis: {e}")
        else:
            if rpm_counts.empty:
                logger.warning("RPM counts not available, skipping enrichment analysis")
            if not input_groups:
                logger.warning("No input groups identified, skipping RPM enrichment analysis")
        
        # Create fixed color palette for groups (for consistent colors across all PCA plots)
        all_groups = meta['group'].unique().tolist()
        group_palette = create_group_color_palette(all_groups)
        logger.info(f"Created color palette for {len(all_groups)} groups")
        
        # Calculate correlation statistics
        df_corr = calculate_correlation_statistics(dds, df_norm)
        
        # Save correlation statistics
        corr_stats_file = output_dir / "correlation_statistics.tsv"
        df_corr.to_csv(corr_stats_file, sep='\t')
        logger.info(f"Saved correlation statistics to {corr_stats_file}")
        
        # Create correlation scatter plots
        corr_plots_file = output_dir / "correlation_scatterPlot.pdf"
        create_correlation_scatter_plots(dds, df_norm, corr_plots_file)
        
        # Perform three PCA analyses
        # 1. PCA with all features
        logger.info("Performing PCA with all features...")
        pca_all, pca_df_all = perform_pca(df_norm, meta, top_n_features=None, 
                                          exclude_input=False, config=config)
        plot_pca(pca_all, pca_df_all, 
                output_dir / "pca_all_features.png",
                title="PCA - All Features",
                group_palette=group_palette, config=config)
        
        # 2. PCA with top variable features
        pca_top_n = config.get("analysis", {}).get("pca", {}).get("top_n_features", 500)
        logger.info(f"Performing PCA with top {pca_top_n} variable features...")
        pca_top500, pca_df_top500 = perform_pca(df_norm, meta, top_n_features=pca_top_n,
                                                exclude_input=False, config=config)
        plot_pca(pca_top500, pca_df_top500,
                output_dir / "pca_top500_features.png",
                title=f"PCA - Top {pca_top_n} Variable Features",
                group_palette=group_palette, config=config)
        
        # 3. PCA without input samples
        logger.info("Performing PCA without input samples...")
        pca_no_input, pca_df_no_input = perform_pca(df_norm, meta, top_n_features=pca_top_n,
                                                    exclude_input=True,
                                                    input_groups=input_groups, config=config)
        plot_pca(pca_no_input, pca_df_no_input,
                output_dir / "pca_no_input_samples.png",
                title=f"PCA - Top {pca_top_n} Features (No Input Samples)",
                group_palette=group_palette, config=config)
        
        # Perform differential expression analysis
        if de_target and de_background:
            logger.info(f"Performing differential expression: {de_target} vs {de_background}")
            try:
                de_results = perform_differential_expression(
                    df_norm, meta, de_target, de_background
                )
                
                # Save differential expression results
                de_file = output_dir / f"DE_{de_target}_vs_{de_background}.tsv"
                de_results.to_csv(de_file, sep='\t')
                logger.info(f"Saved differential expression results to {de_file}")
                
            except Exception as e:
                logger.warning(f"Failed to perform differential expression analysis: {e}")
        else:
            logger.warning("Skipping differential expression analysis: "
                          "Could not determine target and background groups")
        
        logger.info("Statistical analysis completed successfully!")
        return True
        
    except Exception as e:
        logger.error(f"Error during statistical analysis: {e}", exc_info=True)
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
        help='Path to merged.counts.tsv file'
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
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.verbose)
    
    # Run analysis
    success = run_statistical_analysis(
        counts_file=args.counts,
        sample_info_file=args.sample_info,
        output_dir=args.output_dir,
        config_file=args.config,
        n_cpus=args.n_cpus,
        verbose=args.verbose
    )
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()


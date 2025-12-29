#!/usr/bin/env python3
import pandas as pd
import numpy as np
import logging
from pathlib import Path
from Bio.Seq import Seq
from Bio import BiopythonWarning
from typing import Dict, Tuple, Optional, List, Any
import sys
import warnings
import multiprocessing as mp
from functools import partial
import os

def setup_logging(output_dir: Path = None, sample_name: str = None, log_file: Path = None) -> logging.Logger:
    """
    Set up logging configuration. If log_file is provided, use it. Otherwise, use output_dir/sample_name.count.log.
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # Remove all handlers if already set (to avoid duplicate logs)
    if logger.hasHandlers():
        logger.handlers.clear()

    if log_file is None:
        if output_dir is not None and sample_name is not None:
            log_file = output_dir / f"{sample_name}.count.log"
        else:
            log_file = None

    # Create file handler if log_file is provided
    if log_file is not None:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # Always add console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger

def extract_variable_region(seq: str, flank_5p: str, flank_3p: str) -> str:
    """
    Extract variable region between flanking sequences.
    
    Args:
        seq (str): Full sequence
        flank_5p (str): 5' flanking sequence
        flank_3p (str): 3' flanking sequence
        
    Returns:
        str: Extracted variable region or None if flanks not found
    """
    start = seq.find(flank_5p)
    end = seq.find(flank_3p)

    if start != -1 and end != -1:
        start += len(flank_5p)
        return seq[start:end]
    return None

def translate(seq: str) -> str:
    """
    Translate DNA sequence to protein sequence.
    
    Args:
        seq (str): DNA sequence
        
    Returns:
        str: Translated protein sequence or None if translation fails
    """
    try:
        normalized = seq.upper().replace('U', 'T')
        remainder = len(normalized) % 3

        if remainder != 0:
            normalized = normalized[:len(normalized) - remainder]

        if len(normalized) < 3:
            return None

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            return str(Seq(normalized).translate(to_stop=True))
    except Exception:
        return None

def _process_chunk_worker(chunk_lines: List[str], flank_5p: str, flank_3p: str) -> Tuple[List[Dict[str, Any]], Dict[str, int]]:
    """
    Worker function to process a chunk of SAM file lines.
    This function is called by multiprocessing workers.
    
    Args:
        chunk_lines: List of SAM file lines to process
        flank_5p: 5' flanking sequence
        flank_3p: 3' flanking sequence
        
    Returns:
        Tuple of (processed_rows, stats_dict)
    """
    import re
    
    def parse_cigar(cigar: str) -> Dict[str, int]:
        stats = {
            'matches': 0,
            'insertions': 0,
            'deletions': 0,
            'soft_clips': 0,
            'hard_clips': 0,
            'skips': 0,
            'padding': 0
        }
        cigar_parts = re.findall(r'(\d+)([MIDNSHP=XB])', cigar)
        for length, op in cigar_parts:
            length = int(length)
            if op in {'M', '=', 'X'}:
                stats['matches'] += length
            elif op == 'I':
                stats['insertions'] += length
            elif op == 'D':
                stats['deletions'] += length
            elif op == 'S':
                stats['soft_clips'] += length
            elif op == 'H':
                stats['hard_clips'] += length
            elif op == 'N':
                stats['skips'] += length
            elif op == 'P':
                stats['padding'] += length
        return stats
    
    total_reads = 0
    unmapped_reads = 0
    mapped_reads = 0
    reads_with_flanks = 0
    valid_peptides = 0
    processed_rows: List[Dict[str, Any]] = []
    
    for line in chunk_lines:
        if line.startswith("@"):
            continue
        total_reads += 1
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 10:
            continue
        flag = int(parts[1])
        if (flag & 4) != 0:
            unmapped_reads += 1
            continue
        
        mapped_reads += 1
        seq = parts[9]
        variable_seq = extract_variable_region(seq, flank_5p, flank_3p)
        if variable_seq is None:
            continue
        reads_with_flanks += 1
        
        peptide = translate(variable_seq)
        if peptide is None:
            continue
        valid_peptides += 1
        
        cigar_stats = parse_cigar(parts[5])
        processed_rows.append({
            'peptide': peptide,
            'variable_seq': variable_seq,
            'insertions': cigar_stats['insertions'],
            'deletions': cigar_stats['deletions'],
            'matches': cigar_stats['matches']
        })
    
    stats = {
        'total_reads': total_reads,
        'mapped_reads': mapped_reads,
        'unmapped_reads': unmapped_reads,
        'reads_with_flanks': reads_with_flanks,
        'valid_peptides': valid_peptides
    }
    
    return processed_rows, stats

def process_sam_file(sam_file: Path, config: Dict, logger: Optional[logging.Logger] = None, n_cpus: Optional[int] = None) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Process SAM file and extract variable regions using multiprocessing.
    
    Args:
        sam_file (Path): Path to SAM file
        config (Dict): Configuration dictionary containing flanking sequences
        logger (Optional[logging.Logger]): Logger instance
        n_cpus (Optional[int]): Number of CPUs to use. If None, uses all available CPUs.
        
    Returns:
        Tuple[pd.DataFrame, Dict[str, int]]: 
            - DataFrame with processed sequences
            - Dictionary with processing statistics
    """
    # Get flanking sequences from config
    flank_5p = config['flanking_sequences']['flank_5p']
    flank_3p = config['flanking_sequences']['flank_3p']
    
    # Determine number of CPUs to use
    if n_cpus is None:
        n_cpus = max(1, os.cpu_count() or 1)
    else:
        n_cpus = max(1, min(n_cpus, os.cpu_count() or 1))
    
    # Read all non-header lines into memory (chunked for large files)
    if logger:
        logger.info(f"Reading SAM file with {n_cpus} CPU(s)...")
    all_lines = []
    header_lines = []
    
    with open(sam_file, "r") as f:
        for line in f:
            if line.startswith("@"):
                header_lines.append(line)
            else:
                all_lines.append(line)
    
    total_lines = len(all_lines)
    if logger:
        logger.info(f"Total reads to process: {total_lines:,}")
    
    # Determine chunk size (aim for ~100k-500k reads per chunk)
    chunk_size = max(100_000, total_lines // (n_cpus * 4))
    chunks = [all_lines[i:i + chunk_size] for i in range(0, len(all_lines), chunk_size)]
    num_chunks = len(chunks)
    
    if logger:
        logger.info(f"Processing {num_chunks} chunk(s) with ~{chunk_size:,} reads per chunk")
    
    # Create worker function with fixed parameters
    worker_func = partial(_process_chunk_worker, flank_5p=flank_5p, flank_3p=flank_3p)
    
    # Process chunks in parallel
    processed_rows: List[Dict[str, Any]] = []
    total_reads = 0
    unmapped_reads = 0
    mapped_reads = 0
    reads_with_flanks = 0
    valid_peptides = 0
    
    if n_cpus > 1 and num_chunks > 1:
        # Use multiprocessing
        with mp.Pool(processes=n_cpus) as pool:
            results = pool.map(worker_func, chunks)
        
        # Merge results from all chunks
        for chunk_rows, chunk_stats in results:
            processed_rows.extend(chunk_rows)
            total_reads += chunk_stats['total_reads']
            unmapped_reads += chunk_stats['unmapped_reads']
            mapped_reads += chunk_stats['mapped_reads']
            reads_with_flanks += chunk_stats['reads_with_flanks']
            valid_peptides += chunk_stats['valid_peptides']
            
            if logger:
                logger.info(
                    "Processed chunk: %s reads (%s mapped, %s valid peptides)",
                    f"{chunk_stats['total_reads']:,}",
                    f"{chunk_stats['mapped_reads']:,}",
                    f"{chunk_stats['valid_peptides']:,}"
                )
    else:
        # Single-threaded fallback
        if logger:
            logger.info("Using single-threaded processing")
        chunk_rows, chunk_stats = worker_func(chunks[0])
        processed_rows.extend(chunk_rows)
        total_reads = chunk_stats['total_reads']
        unmapped_reads = chunk_stats['unmapped_reads']
        mapped_reads = chunk_stats['mapped_reads']
        reads_with_flanks = chunk_stats['reads_with_flanks']
        valid_peptides = chunk_stats['valid_peptides']
    
    if logger:
        logger.info(f"Total processed: {total_reads:,} reads, {mapped_reads:,} mapped, {valid_peptides:,} valid peptides")
    
    # Create DataFrame from all processed rows
    df = pd.DataFrame(processed_rows)
    
    if df.empty:
        variant_stats = {
            'total_insertions': 0,
            'total_deletions': 0,
            'reads_with_insertions': 0,
            'reads_with_deletions': 0
        }
    else:
        variant_stats = {
            'total_insertions': df['insertions'].sum(),
            'total_deletions': df['deletions'].sum(),
            'reads_with_insertions': (df['insertions'] > 0).sum(),
            'reads_with_deletions': (df['deletions'] > 0).sum()
        }
    
    stats = {
        'total_reads': total_reads,
        'mapped_reads': mapped_reads,
        'unmapped_reads': unmapped_reads,
        'reads_with_flanks': reads_with_flanks,
        'valid_peptides': valid_peptides,
        'mapping_rate': mapped_reads / total_reads if total_reads > 0 else 0,
        'flank_detection_rate': reads_with_flanks / mapped_reads if mapped_reads > 0 else 0,
        'translation_success_rate': valid_peptides / reads_with_flanks if reads_with_flanks > 0 else 0,
        'variant_stats': variant_stats
    }
    
    return df, stats

def merge_with_reference(df: pd.DataFrame, reference_file: Path) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Merge processed reads with reference library.
    
    Args:
        df (pd.DataFrame): DataFrame with processed reads
        reference_file (Path): Path to reference library file
        
    Returns:
        Tuple[pd.DataFrame, Dict[str, int]]: 
            - Merged DataFrame
            - Dictionary with merging statistics
    """
    # Read reference file and keep only essential columns
    ref_df = pd.read_csv(reference_file)[['peptide', 'ID_WLG']]
    
    # Get unique peptides before merging
    unique_peptides = df['peptide'].nunique()
    
    # Merge on peptide sequence
    df_merged = pd.merge(df, ref_df, on='peptide', how='left')
    
    # Calculate statistics
    stats = {
        'unique_peptides': unique_peptides,
        'unique_peptides_in_ref': ref_df['peptide'].nunique(),
        'peptides_in_both': df_merged['ID_WLG'].notna().sum(),
        'peptides_only_in_reads': df_merged['ID_WLG'].isna().sum()
    }
    
    return df_merged, stats

def calculate_counts(df: pd.DataFrame, total_reads: int) -> Tuple[pd.DataFrame, Dict[str, int]]:
    # Calculate peptide abundance (how many times each peptide appears)
    peptide_counts = df['peptide'].value_counts()
    df['count'] = df['peptide'].map(peptide_counts)
    df['RPM'] = df['count'] / total_reads * 1_000_000
    # Add log2(RPM + 1) transformation
    df['log2_RPM'] = np.log2(df['RPM'] + 1)

    # Reorder columns as desired
    df_out = df[['ID_WLG', 'peptide', 'variable_seq', 'count', 'RPM', 'log2_RPM', 'insertions', 'deletions', 'matches']].copy()

    # Stats (optional, can be adjusted)
    assigned = df_out[df_out['ID_WLG'].notna()]['count'].sum()
    unassigned = df_out[df_out['ID_WLG'].isna()]['count'].sum()
    stats = {
        'total': len(df_out),
        'assigned': assigned,
        'unassigned': unassigned,
        'unique_variants': df_out['peptide'].nunique(),
        'max_reads_per_variant': df_out['count'].max(),
        'min_reads_per_variant': df_out['count'].min(),
        'mean_reads_per_variant': df_out['count'].mean()
    }
    return df_out, stats

def main(
    sam_file: Path,
    reference_file: Path,
    config: Dict,
    output_file: Path,
    log_file: Path = None,
    logger: Optional[logging.Logger] = None,
    n_cpus: Optional[int] = None
) -> None:
    """
    Main function to process SAM file and generate counts.
    Accepts an optional log_file argument for unified logging.
    
    Args:
        sam_file: Path to SAM file
        reference_file: Path to reference library file
        config: Configuration dictionary
        output_file: Path to output file
        log_file: Optional log file path
        logger: Optional logger instance
        n_cpus: Number of CPUs to use for parallel processing. If None, uses all available CPUs.
    """
    try:
        # Set up logging
        output_dir = output_file.parent
        sample_name = sam_file.stem.split('.')[0]  # Get sample name from SAM file
        if logger is None:
            logger = setup_logging(output_dir, sample_name, log_file=log_file)
        else:
            logger.info(f"Using existing logger; count output will be appended to {log_file if log_file else 'configured handlers'}.")
        
        logger.info(f"Starting variant counting for sample: {sample_name}")
        logger.info(f"Input SAM file: {sam_file}")
        logger.info(f"Reference library: {reference_file}")
        logger.info(f"Output file: {output_file}")
        logger.info("Step 1/4: Processing SAM file and extracting variable regions...")
        
        # Process SAM file with multiprocessing
        df, sam_stats = process_sam_file(sam_file, config, logger=logger, n_cpus=n_cpus)
        logger.info("Step 1/4 complete.")
        
        # Log SAM processing statistics
        logger.info("SAM File Processing Statistics:")
        logger.info(f"Total reads in SAM file: {sam_stats['total_reads']:,}")
        logger.info(f"Mapped reads: {sam_stats['mapped_reads']:,} ({sam_stats['mapping_rate']:.2%})")
        logger.info(f"Unmapped reads: {sam_stats['unmapped_reads']:,}")
        logger.info(f"Reads with flanking sequences: {sam_stats['reads_with_flanks']:,} ({sam_stats['flank_detection_rate']:.2%})")
        logger.info(f"Valid peptides after translation: {sam_stats['valid_peptides']:,} ({sam_stats['translation_success_rate']:.2%})")
        
        logger.info("Step 2/4: Merging reads with reference library...")
        df_merged, merge_stats = merge_with_reference(df, reference_file)
        logger.info("Step 2/4 complete.")

        # Ensure all unassigned ID_WLG are labeled as 'Unassigned'
        df_merged['ID_WLG'] = df_merged['ID_WLG'].fillna('Unassigned')
        
        # Log merging statistics
        logger.info("\nReference Library Merging Statistics:")
        logger.info(f"Unique peptides in reads: {merge_stats['unique_peptides']:,}")
        logger.info(f"Unique peptides in reference: {merge_stats['unique_peptides_in_ref']:,}")
        logger.info(f"Peptides found in both: {merge_stats['peptides_in_both']:,}")
        logger.info(f"Peptides only in reads: {merge_stats['peptides_only_in_reads']:,}")
        
        logger.info("Step 3/4: Calculating counts and abundance metrics...")
        df_final, count_stats = calculate_counts(df_merged, sam_stats['total_reads'])
        logger.info("Step 3/4 complete.")

        # Calculate correct assigned/unassigned read counts and percentages
        total_sequences = len(df_final)
        assigned_sequences = df_final[df_final['ID_WLG'].notna() & (df_final['ID_WLG'] != 'Unassigned')].shape[0]
        unassigned_sequences = df_final[(df_final['ID_WLG'].isna()) | (df_final['ID_WLG'] == 'Unassigned')].shape[0]

        # Log count statistics with correct percentages
        logger.info("\nVariant Count Statistics:")
        logger.info(f"Total sequences: {total_sequences:,}")
        logger.info(f"Assigned sequences: {assigned_sequences:,} ({assigned_sequences/total_sequences:.2%})")
        logger.info(f"Unassigned sequences: {unassigned_sequences:,} ({unassigned_sequences/total_sequences:.2%})")
        logger.info(f"Unique variants found: {count_stats['unique_variants']:,}")
        logger.info(f"Reads per variant - Max: {count_stats['max_reads_per_variant']:,}, Min: {count_stats['min_reads_per_variant']:,}, Mean: {count_stats['mean_reads_per_variant']:.2f}")
        
        # Prepare final output: drop duplicate peptides, sort by RPM, reset index
        df_final = df_final.drop_duplicates(subset=['peptide'], keep='first')
        df_final = df_final.sort_values(by='RPM', ascending=False)
        df_final = df_final.reset_index(drop=True)
        
        logger.info("Step 4/4: Writing final count table to disk...")
        df_final.to_csv(output_file, index=False)
        logger.info(f"\nResults saved to {output_file}")
        if log_file:
            logger.info(f"Log file saved to {log_file}")
        
    except Exception as e:
        logger.error(f"Error processing SAM file: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Process SAM file and count variants')
    parser.add_argument('sam_file', type=Path, help='Path to SAM file')
    parser.add_argument('reference_file', type=Path, help='Path to reference library file')
    parser.add_argument('config_file', type=Path, help='Path to config file')
    parser.add_argument('--output', type=Path, help='Path to output file')
    
    args = parser.parse_args()
    
    # Read config
    import json
    with open(args.config_file) as f:
        config = json.load(f)
    
    # Set default output path if not provided
    if not args.output:
        args.output = args.sam_file.parent / f"{args.sam_file.stem}.counts.tsv"
    
    main(args.sam_file, args.reference_file, config, args.output)
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
import os
import tempfile
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

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

def process_sam_file(sam_file: Path, config: Dict, n_threads: Optional[int] = None, logger: Optional[logging.Logger] = None) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Process SAM file and extract variable regions using memory-efficient chunked processing.
    
    Uses temporary Parquet files to avoid loading all data into memory at once.
    
    Args:
        sam_file (Path): Path to SAM file
        config (Dict): Configuration dictionary containing flanking sequences
        n_threads (Optional[int]): Not used (kept for API compatibility)
        logger (Optional[logging.Logger]): Logger instance
        
    Returns:
        Tuple[pd.DataFrame, Dict[str, int]]: 
            - DataFrame with processed sequences
            - Dictionary with processing statistics
    """
    # Get flanking sequences from config
    flank_5p = config['flanking_sequences']['flank_5p']
    flank_3p = config['flanking_sequences']['flank_3p']
    
    # Memory-efficient chunked processing parameters
    chunk_size = 500_000  # Process 100k valid reads per chunk
    progress_interval = 1_000_000
    
    # Initialize counters
    total_reads = 0
    unmapped_reads = 0
    mapped_reads = 0
    reads_with_flanks = 0
    valid_peptides = 0
    processed_rows: List[Dict[str, Any]] = []
    chunk_count = 0
    
    # Create temporary directory for Parquet chunks
    temp_dir = Path(tempfile.mkdtemp(prefix="capscreen_chunks_"))
    chunk_files: List[Path] = []
    
    if logger:
        logger.info("Using memory-efficient chunked processing (chunk size: %s reads)", f"{chunk_size:,}")
        logger.info("Temporary chunk directory: %s", temp_dir)
    
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
    
    def write_chunk_to_parquet(rows: List[Dict[str, Any]], chunk_num: int) -> Path:
        """Write a chunk of processed rows to a Parquet file (or CSV if Parquet unavailable)."""
        chunk_df = pd.DataFrame(rows)
        try:
            # Try Parquet first (more efficient)
            chunk_file = temp_dir / f"chunk_{chunk_num:06d}.parquet"
            chunk_df.to_parquet(chunk_file, index=False, compression='snappy')
            return chunk_file
        except (ImportError, AttributeError) as e:
            # Fallback to CSV if Parquet is not available
            if logger:
                if chunk_num == 1:  # Only log once
                    logger.warning(
                        "Parquet support not available (pyarrow not installed), falling back to CSV format. "
                        "This will use more disk space. Consider installing pyarrow for better performance."
                    )
            chunk_file = temp_dir / f"chunk_{chunk_num:06d}.csv.gz"
            chunk_df.to_csv(chunk_file, index=False, compression='gzip')
            return chunk_file
    
    try:
        with open(sam_file, "r") as f:
            for line_number, line in enumerate(f, start=1):
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
                
                # Write chunk to disk when it reaches chunk_size
                if len(processed_rows) >= chunk_size:
                    chunk_count += 1
                    chunk_file = write_chunk_to_parquet(processed_rows, chunk_count)
                    chunk_files.append(chunk_file)
                    if logger:
                        logger.debug(
                            "Wrote chunk %s to disk (%s rows, %s total chunks)",
                            chunk_count,
                            f"{len(processed_rows):,}",
                            len(chunk_files)
                        )
                    processed_rows = []  # Clear memory
                
                if logger and total_reads % progress_interval == 0:
                    logger.info(
                        "Processed %s reads (%s mapped, %s valid peptides, %s chunks written) ...",
                        f"{total_reads:,}",
                        f"{mapped_reads:,}",
                        f"{valid_peptides:,}",
                        len(chunk_files)
                    )
        
        # Write final chunk if there are remaining rows
        if processed_rows:
            chunk_count += 1
            chunk_file = write_chunk_to_parquet(processed_rows, chunk_count)
            chunk_files.append(chunk_file)
            if logger:
                logger.debug(
                    "Wrote final chunk %s to disk (%s rows)",
                    chunk_count,
                    f"{len(processed_rows):,}"
                )
            processed_rows = []  # Clear memory
        
        # Combine all chunks into final DataFrame (sequential reading)
        if logger:
            logger.info("Reading %s chunk file(s) and combining into DataFrame...", len(chunk_files))
        
        if chunk_files:
            df = _combine_chunks_sequential(chunk_files, logger)
            if logger:
                logger.info("Successfully combined chunks into DataFrame with %s rows", f"{len(df):,}")
        else:
            df = pd.DataFrame(columns=['peptide', 'variable_seq', 'insertions', 'deletions', 'matches'])
            if logger:
                logger.warning("No valid reads found in SAM file")
    
    finally:
        # Clean up temporary directory
        if logger:
            logger.debug("Cleaning up temporary chunk directory: %s", temp_dir)
        try:
            shutil.rmtree(temp_dir)
            if logger:
                logger.debug("Temporary directory cleaned up successfully")
        except Exception as e:
            if logger:
                logger.warning("Failed to clean up temporary directory %s: %s", temp_dir, e)
    
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


def _read_chunk_file(chunk_file: Path) -> pd.DataFrame:
    """
    Read a single chunk file (Parquet or CSV).
    
    Args:
        chunk_file: Path to chunk file
        
    Returns:
        DataFrame with chunk data
        
    Raises:
        Exception: If file cannot be read
    """
    try:
        if chunk_file.suffix == '.parquet':
            return pd.read_parquet(chunk_file)
        else:
            return pd.read_csv(chunk_file, compression='gzip')
    except Exception as e:
        raise Exception(f"Error reading chunk file {chunk_file.name}: {e}") from e


def _combine_chunks_parallel(chunk_files: List[Path], n_threads: int, logger: Optional[logging.Logger] = None) -> pd.DataFrame:
    """
    Read and combine chunks in parallel using threading.
    
    This function uses ThreadPoolExecutor to read multiple chunk files concurrently,
    which is effective for I/O-bound operations like reading from disk.
    
    Args:
        chunk_files: List of chunk file paths
        n_threads: Number of threads to use for parallel reading (not currently used)
        logger: Optional logger instance
        
    Returns:
        Combined DataFrame with all chunks
        
    Raises:
        Exception: If any chunk file cannot be read
    """
    if not chunk_files:
        return pd.DataFrame(columns=['peptide', 'variable_seq', 'insertions', 'deletions', 'matches'])
    
    if logger:
        logger.info(f"Reading {len(chunk_files)} chunks in parallel using {n_threads} threads...")
    
    df_parts = []
    errors = []
    
    try:
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            # Submit all read tasks
            future_to_file = {executor.submit(_read_chunk_file, f): f for f in chunk_files}
            
            # Collect results as they complete
            completed = 0
            for future in as_completed(future_to_file):
                chunk_file = future_to_file[future]
                completed += 1
                try:
                    df_chunk = future.result()
                    df_parts.append(df_chunk)
                    if logger and completed % max(1, len(chunk_files) // 10) == 0:
                        logger.debug(
                            "Read %s/%s chunks (%s rows in this chunk)...",
                            completed,
                            len(chunk_files),
                            f"{len(df_chunk):,}"
                        )
                except Exception as e:
                    error_msg = f"Error reading chunk {chunk_file.name}: {e}"
                    errors.append(error_msg)
                    if logger:
                        logger.error(error_msg)
                    # Continue processing other chunks, but track the error
    
    except Exception as e:
        if logger:
            logger.error(f"Fatal error in parallel chunk reading: {e}", exc_info=True)
        raise
    
    # Check if we had any errors
    if errors:
        error_summary = f"Failed to read {len(errors)} chunk(s) out of {len(chunk_files)}"
        if logger:
            logger.error(error_summary)
        raise Exception(f"{error_summary}. First error: {errors[0]}")
    
    # Combine all chunks
    if df_parts:
        if logger:
            logger.debug(f"Combining {len(df_parts)} chunk DataFrames...")
        df = pd.concat(df_parts, ignore_index=True)
        if logger:
            logger.info(f"Successfully combined {len(df_parts)} chunks into DataFrame with {len(df):,} rows")
        return df
    else:
        if logger:
            logger.warning("No chunks were successfully read")
        return pd.DataFrame(columns=['peptide', 'variable_seq', 'insertions', 'deletions', 'matches'])


def _combine_chunks_sequential(chunk_files: List[Path], logger: Optional[logging.Logger] = None) -> pd.DataFrame:
    """
    Read and combine chunks sequentially (original behavior).
    
    This is the fallback method when parallel reading is not requested or not available.
    
    Args:
        chunk_files: List of chunk file paths
        logger: Optional logger instance
        
    Returns:
        Combined DataFrame with all chunks
    """
    if not chunk_files:
        return pd.DataFrame(columns=['peptide', 'variable_seq', 'insertions', 'deletions', 'matches'])
    
    # Determine file format from first chunk
    use_parquet = chunk_files[0].suffix == '.parquet'
    
    # Read chunks sequentially
    if len(chunk_files) <= 50:
        # Small number of chunks: read all at once
        if use_parquet:
            df = pd.concat([pd.read_parquet(f) for f in chunk_files], ignore_index=True)
        else:
            df = pd.concat([pd.read_csv(f, compression='gzip') for f in chunk_files], ignore_index=True)
    else:
        # Large number of chunks: read in batches
        batch_size = 50
        df_parts = []
        for i in range(0, len(chunk_files), batch_size):
            batch = chunk_files[i:i + batch_size]
            if use_parquet:
                batch_df = pd.concat([pd.read_parquet(f) for f in batch], ignore_index=True)
            else:
                batch_df = pd.concat([pd.read_csv(f, compression='gzip') for f in batch], ignore_index=True)
            df_parts.append(batch_df)
            if logger and (i // batch_size + 1) % 10 == 0:
                logger.debug("Read %s/%s chunk batches...", (i // batch_size + 1), (len(chunk_files) // batch_size + 1))
        df = pd.concat(df_parts, ignore_index=True)
    
    return df


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
    n_threads: Optional[int] = None
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
        n_threads: Not used (kept for API compatibility). Threads are read from config['threads'] if available.
    """
    try:
        # Set up logging
        output_dir = output_file.parent
        sample_name = sam_file.stem.split('.')[0]  # Get sample name from SAM file
        if logger is None:
            logger = setup_logging(output_dir, sample_name, log_file=log_file)
        else:
            logger.info(f"Using existing logger; count output will be appended to {log_file if log_file else 'configured handlers'}.")
        
        # Determine number of threads: use n_threads if provided, otherwise read from config
        if n_threads is None:
            n_threads = config.get('threads')
            if n_threads is not None:
                logger.info(f"Using threads from config file: {n_threads}")
        
        logger.info(f"Starting variant counting for sample: {sample_name}")
        logger.info(f"Input SAM file: {sam_file}")
        logger.info(f"Reference library: {reference_file}")
        logger.info(f"Output file: {output_file}")
        logger.info("Step 1/4: Processing SAM file and extracting variable regions...")
        
        # Process SAM file
        df, sam_stats = process_sam_file(sam_file, config, n_threads=n_threads, logger=logger)
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
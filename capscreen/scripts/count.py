#!/usr/bin/env python3
import pandas as pd
import logging
from pathlib import Path
from Bio.Seq import Seq
from typing import Dict, Tuple
import sys

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
        return str(Seq(seq).translate(to_stop=True))
    except:
        return None

def process_sam_file(sam_file: Path, config: Dict) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Process SAM file and extract variable regions.
    
    Args:
        sam_file (Path): Path to SAM file
        config (Dict): Configuration dictionary containing flanking sequences
        
    Returns:
        Tuple[pd.DataFrame, Dict[str, int]]: 
            - DataFrame with processed sequences
            - Dictionary with processing statistics
    """
    # Get flanking sequences from config
    flank_5p = config['flanking_sequences']['flank_5p']
    flank_3p = config['flanking_sequences']['flank_3p']
    
    # Read SAM file and extract mapped reads
    reads = []
    total_reads = 0
    unmapped_reads = 0
    
    with open(sam_file, "r") as f:
        for line in f:
            if not line.startswith("@"):
                total_reads += 1
                parts = line.strip().split("\t")
                flag = int(parts[1])
                if (flag & 4) == 0:  # only mapped reads
                    seq = parts[9]
                    cigar = parts[5]
                    reads.append({
                        'full_seq': seq,
                        'cigar': cigar,
                        'reference': parts[2],  # Reference sequence name
                        'position': int(parts[3]),  # 1-based position
                        'mapping_quality': int(parts[4])  # Mapping quality
                    })
                else:
                    unmapped_reads += 1
    
    # Create DataFrame
    df = pd.DataFrame(reads)
    
    # Extract variable regions
    df['variable_seq'] = df['full_seq'].apply(lambda s: extract_variable_region(s, flank_5p, flank_3p))
    reads_with_flanks = df['variable_seq'].notna().sum()
    df = df.dropna(subset=['variable_seq'])  # remove reads where flanks not found
    
    # Translate to peptides
    df['peptide'] = df['variable_seq'].apply(translate)
    valid_peptides = df['peptide'].notna().sum()
    df = df.dropna(subset=['peptide'])  # remove frameshifted or invalid
    
    # Parse CIGAR string to get variant information
    def parse_cigar(cigar: str) -> Dict:
        import re
        stats = {
            'matches': 0,
            'insertions': 0,
            'deletions': 0,
            'soft_clips': 0,
            'hard_clips': 0,
            'skips': 0,
            'padding': 0
        }
        
        # Parse CIGAR string
        cigar_parts = re.findall(r'(\d+)([MIDNSHP=XB])', cigar)
        for length, op in cigar_parts:
            length = int(length)
            if op == 'M' or op == '=' or op == 'X':
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
    
    # Add variant information to DataFrame
    df['cigar_stats'] = df['cigar'].apply(parse_cigar)
    df['insertions'] = df['cigar_stats'].apply(lambda x: x['insertions'])
    df['deletions'] = df['cigar_stats'].apply(lambda x: x['deletions'])
    df['matches'] = df['cigar_stats'].apply(lambda x: x['matches'])
    
    # Compile statistics
    stats = {
        'total_reads': total_reads,
        'mapped_reads': len(reads),
        'unmapped_reads': unmapped_reads,
        'reads_with_flanks': reads_with_flanks,
        'valid_peptides': valid_peptides,
        'mapping_rate': len(reads) / total_reads if total_reads > 0 else 0,
        'flank_detection_rate': reads_with_flanks / len(reads) if len(reads) > 0 else 0,
        'translation_success_rate': valid_peptides / reads_with_flanks if reads_with_flanks > 0 else 0,
        'variant_stats': {
            'total_insertions': df['insertions'].sum(),
            'total_deletions': df['deletions'].sum(),
            'reads_with_insertions': (df['insertions'] > 0).sum(),
            'reads_with_deletions': (df['deletions'] > 0).sum()
        }
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

    # Reorder columns as desired
    df_out = df[['ID_WLG', 'peptide', 'variable_seq', 'count', 'RPM', 'insertions', 'deletions', 'matches']].copy()

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

def main(sam_file: Path, reference_file: Path, config: Dict, output_file: Path, log_file: Path = None) -> None:
    """
    Main function to process SAM file and generate counts.
    Accepts an optional log_file argument for unified logging.
    """
    try:
        # Set up logging
        output_dir = output_file.parent
        sample_name = sam_file.stem.split('.')[0]  # Get sample name from SAM file
        logger = setup_logging(output_dir, sample_name, log_file=log_file)
        
        logger.info(f"Starting variant counting for sample: {sample_name}")
        logger.info(f"Input SAM file: {sam_file}")
        logger.info(f"Reference library: {reference_file}")
        logger.info(f"Output file: {output_file}")
        
        # Process SAM file
        df, sam_stats = process_sam_file(sam_file, config)
        
        # Log SAM processing statistics
        logger.info("SAM File Processing Statistics:")
        logger.info(f"Total reads in SAM file: {sam_stats['total_reads']:,}")
        logger.info(f"Mapped reads: {sam_stats['mapped_reads']:,} ({sam_stats['mapping_rate']:.2%})")
        logger.info(f"Unmapped reads: {sam_stats['unmapped_reads']:,}")
        logger.info(f"Reads with flanking sequences: {sam_stats['reads_with_flanks']:,} ({sam_stats['flank_detection_rate']:.2%})")
        logger.info(f"Valid peptides after translation: {sam_stats['valid_peptides']:,} ({sam_stats['translation_success_rate']:.2%})")
        
        # Merge with reference
        df_merged, merge_stats = merge_with_reference(df, reference_file)

        # Ensure all unassigned ID_WLG are labeled as 'Unassigned'
        df_merged['ID_WLG'] = df_merged['ID_WLG'].fillna('Unassigned')
        
        # Calculate variant statistics for unassigned reads only
        unassigned_df = df_merged[df_merged['ID_WLG'].isna()]
        unassigned_variant_stats = {
            'total_insertions': unassigned_df['insertions'].sum(),
            'total_deletions': unassigned_df['deletions'].sum(),
            'reads_with_insertions': (unassigned_df['insertions'] > 0).sum(),
            'reads_with_deletions': (unassigned_df['deletions'] > 0).sum()
        }
        
        # Log variant statistics for unassigned reads
        logger.info("\nVariant Statistics for Unassigned Reads:")
        logger.info(f"Total insertions: {unassigned_variant_stats['total_insertions']:,}")
        logger.info(f"Total deletions: {unassigned_variant_stats['total_deletions']:,}")
        logger.info(f"Reads with insertions: {unassigned_variant_stats['reads_with_insertions']:,}")
        logger.info(f"Reads with deletions: {unassigned_variant_stats['reads_with_deletions']:,}")
        
        # Log merging statistics
        logger.info("\nReference Library Merging Statistics:")
        logger.info(f"Unique peptides in reads: {merge_stats['unique_peptides']:,}")
        logger.info(f"Unique peptides in reference: {merge_stats['unique_peptides_in_ref']:,}")
        logger.info(f"Peptides found in both: {merge_stats['peptides_in_both']:,}")
        logger.info(f"Peptides only in reads: {merge_stats['peptides_only_in_reads']:,}")
        
        # Calculate counts
        df_final, count_stats = calculate_counts(df_merged, sam_stats['total_reads'])

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
        
        # Save results
        df_final.to_csv(output_file)
        logger.info(f"\nResults saved to {output_file}")
        logger.info(f"Log file saved to {output_dir}/{sample_name}.count.log")
        
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
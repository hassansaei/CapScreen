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

DEFAULT_MAX_MISMATCHES = 2
DEFAULT_PRIMER_N = 1
DEFAULT_MAX_VARIABLE_NT = 120


def _flank_params(config: Dict) -> Tuple[str, str, int, int, int]:
    """Return flank sequences and matching parameters from config."""
    fs = config['flanking_sequences']
    max_mm = int(fs.get('max_mismatches', DEFAULT_MAX_MISMATCHES))
    primer_n = int(fs.get('primer_n', DEFAULT_PRIMER_N))
    max_variable_nt = int(fs.get('max_variable_nt', DEFAULT_MAX_VARIABLE_NT))
    if max_mm < 0:
        max_mm = 0
    if primer_n < 0:
        primer_n = 0
    if max_variable_nt < 1:
        max_variable_nt = DEFAULT_MAX_VARIABLE_NT
    return fs['flank_5p'], fs['flank_3p'], max_mm, primer_n, max_variable_nt


def _hamming_le(a: str, b: str, max_mm: int) -> Optional[int]:
    """Return Hamming distance if it is at most max_mm, else None."""
    if len(a) != len(b):
        return None
    mm = 0
    for x, y in zip(a, b):
        if x != y:
            mm += 1
            if mm > max_mm:
                return None
    return mm


def _find_flank(
    seq: str,
    pattern: str,
    max_mm: int,
    start: int = 0,
    end: Optional[int] = None,
    rightmost: bool = False
) -> Optional[Tuple[int, int]]:
    """
    Find pattern in seq allowing up to max_mm substitutions.

    Prefers an exact match. Otherwise returns the match with the fewest
    mismatches (leftmost, or rightmost if requested).
    """
    n = len(pattern)
    if n == 0 or len(seq) < n:
        return None
    last = (len(seq) if end is None else end) - n
    if last < start:
        return None

    best = None
    positions = range(last, start - 1, -1) if rightmost else range(start, last + 1)
    for i in positions:
        mm = _hamming_le(seq[i:i + n], pattern, max_mm)
        if mm is None:
            continue
        if mm == 0:
            return (i, 0)
        if best is None or mm < best[1]:
            best = (i, mm)
    return best


def _strip_leftover_flanks(
    seq: str,
    flank_5p: str,
    flank_3p: str,
    max_mm: int,
    primer_n: int
) -> Tuple[str, Dict[str, Any]]:
    """
    Remove untrimmed 5'/3' copies of the flanks from an extracted insert.

    Allows up to primer_n random bases (Illumina primer N spacer) outside
    each leftover flank. The 3' search looks near the end of the sequence
    so extra bases after a mismatched 3' flank are removed with it.
    """
    info = {
        'leftover_5p_mm': None,
        'leftover_3p_mm': None,
        'primer_n_5p': False,
        'primer_n_3p': False
    }
    if not seq:
        return seq, info

    n5 = len(flank_5p)
    if n5 and len(seq) >= n5:
        for skip in range(min(primer_n, len(seq) - n5) + 1):
            mm = _hamming_le(seq[skip:skip + n5], flank_5p, max_mm)
            if mm is None:
                continue
            info['leftover_5p_mm'] = mm
            info['primer_n_5p'] = skip > 0
            seq = seq[skip + n5:]
            break

    n3 = len(flank_3p)
    if n3 and len(seq) >= n3:
        extra_tail = primer_n + 16
        window_start = max(0, len(seq) - n3 - extra_tail)
        best = None
        for i in range(window_start, len(seq) - n3 + 1):
            mm = _hamming_le(seq[i:i + n3], flank_3p, max_mm)
            if mm is None:
                continue
            if best is None or mm < best[1] or (mm == best[1] and i < best[0]):
                best = (i, mm)
        if best is not None:
            pos, mm = best
            info['leftover_3p_mm'] = mm
            info['primer_n_3p'] = (len(seq) - (pos + n3)) > 0 and (len(seq) - (pos + n3)) <= primer_n
            seq = seq[:pos]
    return seq, info


def extract_variable_region(
    seq: str,
    flank_5p: str,
    flank_3p: str,
    max_mismatches: int = DEFAULT_MAX_MISMATCHES,
    primer_n: int = DEFAULT_PRIMER_N,
    max_variable_nt: int = DEFAULT_MAX_VARIABLE_NT
) -> Tuple[Optional[str], Dict[str, Any]]:
    """
    Extract the variable region between 5' and 3' flanks.

    Flanks may differ from the configured sequences by up to max_mismatches
    substitutions. After the outer flanks are removed, leftover inner copies
    (untrimmed mismatched flanks) and primer N spacers at both ends are
    stripped when the insert is longer than max_variable_nt.
    """
    info = {
        'mm_5p': None,
        'mm_3p': None,
        'leftover_5p_mm': None,
        'leftover_3p_mm': None,
        'primer_n_5p': False,
        'primer_n_3p': False,
        'corrected': False
    }
    if not seq or not flank_5p or not flank_3p:
        return None, info

    seq = seq.upper()
    flank_5p = flank_5p.upper()
    flank_3p = flank_3p.upper()

    start = seq.find(flank_5p)
    end = seq.rfind(flank_3p)
    if start != -1 and end != -1 and end >= start + len(flank_5p):
        interior = seq[start + len(flank_5p):end]
        info['mm_5p'] = 0
        info['mm_3p'] = 0
    else:
        hit5 = _find_flank(seq, flank_5p, max_mismatches, start=0, rightmost=False)
        if hit5 is None:
            return None, info
        start, mm5 = hit5
        interior_start = start + len(flank_5p)
        hit3 = _find_flank(
            seq, flank_3p, max_mismatches, start=interior_start, rightmost=True
        )
        if hit3 is None:
            return None, info
        end, mm3 = hit3
        if end < interior_start:
            return None, info
        interior = seq[interior_start:end]
        info['mm_5p'] = mm5
        info['mm_3p'] = mm3

    if len(interior) > max_variable_nt:
        interior, leftover = _strip_leftover_flanks(
            interior, flank_5p, flank_3p, max_mismatches, primer_n
        )
        info['leftover_5p_mm'] = leftover['leftover_5p_mm']
        info['leftover_3p_mm'] = leftover['leftover_3p_mm']
        info['primer_n_5p'] = leftover['primer_n_5p']
        info['primer_n_3p'] = leftover['primer_n_3p']
        info['corrected'] = (
            leftover['leftover_5p_mm'] is not None
            or leftover['leftover_3p_mm'] is not None
            or leftover['primer_n_5p']
            or leftover['primer_n_3p']
        )

    if info['mm_5p'] or info['mm_3p']:
        info['corrected'] = True

    return interior, info


def _new_flank_counters(max_mm: int) -> Dict[str, Any]:
    size = max_mm + 1
    return {
        'outer_5p': [0] * size,
        'outer_3p': [0] * size,
        'outer_5p_not_found': 0,
        'leftover_5p': [0] * size,
        'leftover_3p': [0] * size,
        'primer_n_5p': 0,
        'primer_n_3p': 0,
        'corrected': 0
    }


def _record_flank_info(counters: Dict[str, Any], info: Dict[str, Any], found: bool) -> None:
    if not found:
        counters['outer_5p_not_found'] += 1
        return
    mm5 = info.get('mm_5p')
    mm3 = info.get('mm_3p')
    if mm5 is not None and mm5 < len(counters['outer_5p']):
        counters['outer_5p'][mm5] += 1
    if mm3 is not None and mm3 < len(counters['outer_3p']):
        counters['outer_3p'][mm3] += 1
    leftover5 = info.get('leftover_5p_mm')
    leftover3 = info.get('leftover_3p_mm')
    if leftover5 is not None and leftover5 < len(counters['leftover_5p']):
        counters['leftover_5p'][leftover5] += 1
    if leftover3 is not None and leftover3 < len(counters['leftover_3p']):
        counters['leftover_3p'][leftover3] += 1
    if info.get('primer_n_5p'):
        counters['primer_n_5p'] += 1
    if info.get('primer_n_3p'):
        counters['primer_n_3p'] += 1
    if info.get('corrected'):
        counters['corrected'] += 1


def _format_mm_hist(counts: List[int]) -> str:
    parts = []
    for i, n in enumerate(counts):
        if not n:
            continue
        label = "mismatch" if i == 1 else "mismatches"
        parts.append("{0} {1}: {2:,}".format(i, label, n))
    return '; '.join(parts) if parts else 'none'


def _log_flank_stats(
    logger: logging.Logger,
    counters: Dict[str, Any],
    max_mm: int,
    primer_n: int,
    mapped_reads: int
) -> None:
    logger.info(
        "Flank matching settings: max_mismatches=%s, primer_n=%s (applied to 5' and 3' ends)",
        max_mm,
        primer_n
    )
    logger.info("Outer 5' flank mismatches: %s; not found=%s", _format_mm_hist(counters['outer_5p']), f"{counters['outer_5p_not_found']:,}")
    logger.info("Outer 3' flank mismatches: %s", _format_mm_hist(counters['outer_3p']))

    n5 = sum(counters['leftover_5p'])
    n3 = sum(counters['leftover_3p'])
    logger.info(
        "Leftover 5' flank corrected: %s reads (%s)",
        f"{n5:,}",
        _format_mm_hist(counters['leftover_5p'])
    )
    logger.info(
        "Leftover 3' flank corrected: %s reads (%s)",
        f"{n3:,}",
        _format_mm_hist(counters['leftover_3p'])
    )
    logger.info(
        "Primer spacer N trimmed: 5' %s reads; 3' %s reads",
        f"{counters['primer_n_5p']:,}",
        f"{counters['primer_n_3p']:,}"
    )
    pct = (counters['corrected'] / mapped_reads * 100) if mapped_reads else 0.0
    logger.info(
        "Reads with flank mismatches found and corrected: %s / %s (%.2f%%)",
        f"{counters['corrected']:,}",
        f"{mapped_reads:,}",
        pct
    )


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
    flank_5p, flank_3p, max_mm, primer_n, max_variable_nt = _flank_params(config)
    flank_counters = _new_flank_counters(max_mm)
    if logger:
        logger.info(
            "Variable-region extraction: 5' and 3' flanks allow up to %s mismatches; "
            "primer_n=%s (Illumina R1/R2 N spacer); leftover flanks trimmed above %s nt",
            max_mm,
            primer_n,
            max_variable_nt
        )

    # Memory-efficient chunked processing parameters
    chunk_size = 500_000
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
                # Skip empty lines
                if not line.strip():
                    continue
                # Validate SAM line structure before counting
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10:
                    continue
                # Count only valid SAM lines
                total_reads += 1
                flag = int(parts[1])
                if (flag & 4) != 0:
                    unmapped_reads += 1
                    continue
                
                mapped_reads += 1
                seq = parts[9]
                variable_seq, flank_info = extract_variable_region(
                    seq,
                    flank_5p,
                    flank_3p,
                    max_mismatches=max_mm,
                    primer_n=primer_n,
                    max_variable_nt=max_variable_nt
                )
                if variable_seq is None:
                    _record_flank_info(flank_counters, flank_info, found=False)
                    continue
                _record_flank_info(flank_counters, flank_info, found=True)
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
        'variant_stats': variant_stats,
        'flank_counters': flank_counters,
        'max_mismatches': max_mm,
        'primer_n': primer_n
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
    # Read reference file and keep only essential columns.
    # Deduplicate by peptide so that a many-to-many merge does not inflate
    # downstream counts when the reference maps multiple ID_WLG to one peptide.
    ref_df = pd.read_csv(reference_file)[['peptide', 'ID_WLG']].drop_duplicates(subset=['peptide'], keep='first')
    
    # Get unique peptides before merging
    unique_peptides = df['peptide'].nunique()
    
    # Merge on peptide sequence (now guaranteed one-to-many: one ref row per peptide)
    df_merged = pd.merge(df, ref_df, on='peptide', how='left')
    
    # Calculate statistics
    unique_peptides_in_ref = ref_df['peptide'].nunique()
    # Get unique peptides that were found in both (peptides that matched reference)
    unique_peptides_found_in_both = df_merged[df_merged['ID_WLG'].notna()]['peptide'].nunique()
    # Calculate peptides in reference but not in input
    peptides_in_ref_not_in_input = unique_peptides_in_ref - unique_peptides_found_in_both
    peptides_in_ref_not_in_input_pct = (peptides_in_ref_not_in_input / unique_peptides_in_ref * 100) if unique_peptides_in_ref > 0 else 0.0
    
    stats = {
        'unique_peptides_in_ref': unique_peptides_in_ref,
        'peptides_in_ref_not_in_input': peptides_in_ref_not_in_input,
        'peptides_in_ref_not_in_input_pct': peptides_in_ref_not_in_input_pct
    }
    
    return df_merged, stats

def calculate_counts(df: pd.DataFrame, mapped_reads: int) -> Tuple[pd.DataFrame, Dict[str, int]]:
    # Count reads per unique (peptide, variable_seq) pair so that different DNA
    # variants translating to the same peptide each get their own correct count.
    df['count'] = df.groupby(['peptide', 'variable_seq'])['peptide'].transform('size')
    df['RPM'] = df['count'] / mapped_reads * 1_000_000
    df['log2_RPM'] = np.log2(df['RPM'] + 1)

    # Reorder columns as desired
    df_out = df[['ID_WLG', 'peptide', 'variable_seq', 'count', 'RPM', 'log2_RPM', 'insertions', 'deletions', 'matches']].copy()

    # Stats
    assigned = df_out[df_out['ID_WLG'].notna()]['count'].sum()
    unassigned = df_out[df_out['ID_WLG'].isna()]['count'].sum()
    stats = {
        'total': len(df_out),
        'assigned': assigned,
        'unassigned': unassigned,
        'unique_variants': df_out[['peptide', 'variable_seq']].drop_duplicates().shape[0],
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
        _log_flank_stats(
            logger,
            sam_stats['flank_counters'],
            sam_stats['max_mismatches'],
            sam_stats['primer_n'],
            sam_stats['mapped_reads']
        )
        
        logger.info("Step 2/4: Merging reads with reference library...")
        df_merged, merge_stats = merge_with_reference(df, reference_file)
        logger.info("Step 2/4 complete.")

        # Ensure all unassigned ID_WLG are labeled as 'Unassigned'
        df_merged['ID_WLG'] = df_merged['ID_WLG'].fillna('Unassigned')
        
        # Log merging statistics
        logger.info("\nReference Library Merging Statistics:")
        logger.info(f"Unique peptides in reference: {merge_stats['unique_peptides_in_ref']:,}")
        logger.info(f"Peptides in reference but not in input: {merge_stats['peptides_in_ref_not_in_input']:,} ({merge_stats['peptides_in_ref_not_in_input_pct']:.2f}%)")
        
        logger.info("Step 3/4: Calculating counts and abundance metrics...")
        # Use mapped_reads for RPM normalization (not total_reads including unmapped)
        df_final, count_stats = calculate_counts(df_merged, sam_stats['mapped_reads'])
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
        
        # Deduplicate by (peptide, variable_seq) so each unique DNA variant is one row
        df_final = df_final.drop_duplicates(subset=['peptide', 'variable_seq'], keep='first')

        logger.info("Step 4/4: Writing final count tables to disk...")
        mapped_reads = sam_stats['mapped_reads']

        # Split assigned and unassigned
        df_assigned = df_final[df_final['ID_WLG'] != 'Unassigned'].copy()
        df_unassigned = df_final[df_final['ID_WLG'] == 'Unassigned'].copy()

        # Aggregate assigned by (ID_WLG, peptide): sum counts, keep the
        # variable_seq from the variant with the highest read count.
        if not df_assigned.empty:
            df_assigned = df_assigned.sort_values('count', ascending=False)
            df_assigned = df_assigned.groupby(['ID_WLG', 'peptide'], as_index=False).agg({
                'variable_seq': 'first',
                'count': 'sum',
                'insertions': 'first',
                'deletions': 'first',
                'matches': 'first'
            })
            df_assigned['RPM'] = df_assigned['count'] / mapped_reads * 1_000_000
            df_assigned['log2_RPM'] = np.log2(df_assigned['RPM'] + 1)
            df_assigned = df_assigned[['ID_WLG', 'peptide', 'variable_seq', 'count', 'RPM', 'log2_RPM', 'insertions', 'deletions', 'matches']]
            df_assigned = df_assigned.sort_values('RPM', ascending=False).reset_index(drop=True)

        # Write assigned counts (main output)
        df_assigned.to_csv(output_file, index=False)
        logger.info(f"Assigned counts ({len(df_assigned):,} variants) saved to {output_file}")

        # Write unassigned counts to a separate file
        unassigned_file = output_file.parent / f"{sample_name}.unassigned.counts.csv"
        if not df_unassigned.empty:
            df_unassigned = df_unassigned.sort_values('RPM', ascending=False).reset_index(drop=True)
            df_unassigned.to_csv(unassigned_file, index=False)
            logger.info(f"Unassigned counts ({len(df_unassigned):,} variants) saved to {unassigned_file}")
        else:
            logger.info("No unassigned peptides found.")

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
        args.output = args.sam_file.parent / f"{args.sam_file.stem}.counts.csv"
    
    main(args.sam_file, args.reference_file, config, args.output)
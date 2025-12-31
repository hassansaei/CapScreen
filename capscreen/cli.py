#!/usr/bin/env python3
import argparse
import csv
import json
import logging
import os
import sys
import shutil
import subprocess
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Any, Tuple, List, Set
import importlib.resources as pkg_resources
from capscreen.version import __version__
from capscreen.scripts import count as count_module
from capscreen.scripts import generate_report as generate_report_module
from capscreen.scripts import alignment as alignment_module
from capscreen.scripts import qc as qc_module
# stat_module imported conditionally - may use separate conda environment
try:
    from capscreen.scripts import stat as stat_module
    STAT_MODULE_AVAILABLE = True
except ImportError:
    STAT_MODULE_AVAILABLE = False
    stat_module = None


def get_stat_python_executable() -> Optional[Path]:
    """
    Find the Python executable for the capscreen-stat conda environment.
    
    Returns:
        Path to Python executable in stat environment, or None if not found
    """
    # Check common conda environment locations
    possible_paths = [
        Path("/opt/conda/envs/capscreen-stat/bin/python"),  # Docker
        Path.home() / ".conda" / "envs" / "capscreen-stat" / "bin" / "python",  # Local conda
        Path(os.environ.get("CONDA_PREFIX", "")) / "../capscreen-stat/bin/python",  # Relative to current env
    ]
    
    # Also check if CONDA_ENVS_PATH is set
    if "CONDA_ENVS_PATH" in os.environ:
        possible_paths.append(Path(os.environ["CONDA_ENVS_PATH"]) / "capscreen-stat" / "bin" / "python")
    
    for path in possible_paths:
        if path.exists() and path.is_file():
            return path
    
    return None


def run_statistical_analysis_with_stat_env(
    counts_file: Path,
    sample_info_file: Path,
    output_dir: Path,
    config_file: Optional[Path] = None,
    n_cpus: Optional[int] = None,
    verbose: bool = False,
    logger_instance: Optional[logging.Logger] = None
) -> bool:
    """
    Run statistical analysis using the capscreen-stat conda environment.
    
    This function uses subprocess to call stat.py with the stat environment's Python,
    which has Python 3.9+ and all required dependencies (pydeseq2, etc.).
    
    Args:
        counts_file: Path to merged.counts.tsv file
        sample_info_file: Path to Sample_info.csv file
        output_dir: Output directory for results
        config_file: Path to config.json file (optional)
        n_cpus: Number of CPUs to use (optional)
        verbose: Enable verbose logging
        logger_instance: Optional logger instance
        
    Returns:
        True if analysis completed successfully, False otherwise
    """
    stat_python = get_stat_python_executable()
    
    if stat_python is None:
        # Fall back to trying direct import (might work if in stat env already)
        if STAT_MODULE_AVAILABLE and stat_module:
            if logger_instance:
                logger_instance.warning("Stat environment not found, trying direct import...")
            try:
                return stat_module.run_statistical_analysis(
                    counts_file=counts_file,
                    sample_info_file=sample_info_file,
                    output_dir=output_dir,
                    config_file=config_file,
                    n_cpus=n_cpus,
                    verbose=verbose,
                    logger_instance=logger_instance
                )
            except Exception as e:
                if logger_instance:
                    logger_instance.error(f"Direct import failed: {e}")
                return False
        else:
            if logger_instance:
                logger_instance.error("Cannot run statistical analysis: capscreen-stat environment not found and stat module not available")
            return False
    
    # Build command to run stat.py as a module using the stat environment's Python
    # First try to find the file path via import (most reliable)
    stat_script_path = None
    try:
        import capscreen.scripts.stat as stat_module_file
        if hasattr(stat_module_file, '__file__'):
            stat_script_path = Path(stat_module_file.__file__)
    except Exception:
        pass
    
    # Use module execution if we can't find the file, otherwise use direct file execution
    if stat_script_path and stat_script_path.exists():
        # Use direct file execution (more explicit)
        cmd = [
            str(stat_python),
            str(stat_script_path),
            "--counts", str(counts_file),
            "--sample-info", str(sample_info_file),
            "--output-dir", str(output_dir),
        ]
    else:
        # Fallback to module execution
        cmd = [
            str(stat_python),
            "-m", "capscreen.scripts.stat",
            "--counts", str(counts_file),
            "--sample-info", str(sample_info_file),
            "--output-dir", str(output_dir),
        ]
    
    if config_file:
        cmd.extend(["--config", str(config_file)])
    
    if n_cpus:
        cmd.extend(["--n-cpus", str(n_cpus)])
    
    if verbose:
        cmd.append("--verbose")
    
    # Run the command
    try:
        if logger_instance:
            logger_instance.info(f"Running statistical analysis using stat environment: {stat_python}")
            logger_instance.debug(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False
        )
        
        # Log output
        if logger_instance:
            if result.stdout:
                logger_instance.info(result.stdout)
            if result.stderr:
                if result.returncode != 0:
                    logger_instance.error(result.stderr)
                else:
                    logger_instance.debug(result.stderr)
        
        return result.returncode == 0
        
    except Exception as e:
        if logger_instance:
            logger_instance.error(f"Failed to run statistical analysis: {e}", exc_info=True)
        return False

# Global Logger
LOG_FORMAT = "%(asctime)s [%(levelname)s] [%(name)s:%(lineno)d] %(message)s"
logger = logging.getLogger("FastQProcessor")

SAMPLE_INFO_REQUIRED_COLUMNS = {
    "sample_id",
    "run_mode",
    "path_FASTQ",
    "group",
    "tech_rep",
    "bio_rep",
    "batch_seq"
}

@dataclass
class SampleInfo:
    sample_id: str
    run_mode: str
    fastq1: Path
    fastq2: Optional[Path]
    group: str
    tech_rep: int
    bio_rep: int
    batch_seq: int

@contextmanager
def sample_file_log_handler(log_path: Optional[Path]):
    """
    Temporarily attach a file handler to the global logger so that
    messages are duplicated into the provided sample-specific log file.
    """
    handler = None
    if log_path:
        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            handler = logging.FileHandler(log_path, mode='w')
            handler.setFormatter(logging.Formatter(LOG_FORMAT))
            logger.addHandler(handler)
        except Exception as exc:
            logger.error(f"Failed to attach sample log handler for {log_path}: {exc}", exc_info=True)
            handler = None
    try:
        yield
    finally:
        if handler:
            logger.removeHandler(handler)
            handler.close()

def detach_log_file_handler(logger_instance: logging.Logger, handler: Optional[logging.Handler]) -> None:
    """
    Remove and close a logging file handler if it exists.
    """
    if handler:
        handler.flush()
        logger_instance.removeHandler(handler)
        handler.close()

def append_sample_log_to_batch(
    batch_log_path: Optional[Path],
    sample_log_path: Path,
    sample_name: str,
    fallback_message: Optional[str] = None,
    copy_sample_log: bool = True
) -> None:
    """
    Append the per-sample pipeline log into the batch log so the combined file
    mirrors the sample logs sequentially.
    """
    if not batch_log_path:
        return
    try:
        if not sample_log_path.exists() and not fallback_message:
            return
        batch_log_path.parent.mkdir(parents=True, exist_ok=True)
        with batch_log_path.open("a") as batch_log:
            batch_log.write(f"\n===== Sample {sample_name} =====\n")
            if copy_sample_log and sample_log_path.exists():
                with sample_log_path.open("r") as sample_log:
                    shutil.copyfileobj(sample_log, batch_log)
                batch_log.write("\n")
            elif fallback_message:
                batch_log.write(f"{fallback_message.strip()}\n")
            else:
                batch_log.write("No log output was generated for this sample.\n")
    except Exception as exc:
        logger.error(
            f"Failed to append log for sample {sample_name} to batch log {batch_log_path}: {exc}",
            exc_info=True
        )

def append_batch_summary(batch_log_path: Optional[Path], lines: List[str], header: str = "Batch Summary") -> None:
    """
    Append a short summary block to the batch log file.
    """
    if not batch_log_path:
        return
    try:
        batch_log_path.parent.mkdir(parents=True, exist_ok=True)
        with batch_log_path.open("a") as batch_log:
            batch_log.write(f"\n===== {header} =====\n")
            for line in lines:
                batch_log.write(f"{line.rstrip()}\n")
    except Exception as exc:
        logger.error(f"Failed to append batch summary to {batch_log_path}: {exc}", exc_info=True)

def determine_keep_intermediate(cli_flag: bool, config: Dict[str, Any]) -> bool:
    """
    Resolve whether intermediate files should be preserved by combining
    the CLI flag with the configuration default.
    """
    if cli_flag:
        return True
    return config.get('keep_intermediate', False)

def infer_fastq2_path(fastq1: Path) -> Path:
    """
    Infer the R2 FASTQ path from the R1 FASTQ path.
    
    Supports multiple naming conventions:
    - Illumina standard: _1.fq.gz / _2.fq.gz or _1.fastq.gz / _2.fastq.gz
    - Common patterns: _R1_ / _R2_, _R1. / _R2., _R1 / _R2
    - Alternative: R1_ / R2_, R1. / R2.
    """
    name = fastq1.name
    replacements = [
        ("_1.fq.gz", "_2.fq.gz"),      # Illumina standard: _1.fq.gz → _2.fq.gz
        ("_1.fastq.gz", "_2.fastq.gz"), # Illumina standard: _1.fastq.gz → _2.fastq.gz
        ("_1.fq", "_2.fq"),             # Without compression: _1.fq → _2.fq
        ("_1.fastq", "_2.fastq"),       # Without compression: _1.fastq → _2.fastq
        ("_1.", "_2."),                 # Generic: _1. → _2. (catches other extensions)
        ("_R1_", "_R2_"),               # Standard: _R1_ → _R2_
        ("_R1.", "_R2."),               # Standard: _R1. → _R2.
        ("_R1", "_R2"),                 # Standard: _R1 → _R2
        ("R1_", "R2_"),                 # Alternative: R1_ → R2_
        ("R1.", "R2.")                  # Alternative: R1. → R2.
    ]
    for old, new in replacements:
        if old in name:
            return fastq1.with_name(name.replace(old, new, 1))
    raise ValueError(f"Unable to infer R2 FASTQ file name from {fastq1}")

def _parse_int_field(value: Optional[str], field_name: str) -> int:
    """
    Parse an integer field from the CSV and raise a ValueError with context on failure.
    """
    if value is None:
        raise ValueError(f"Missing value for {field_name}")
    text = str(value).strip()
    if not text:
        raise ValueError(f"Empty value for {field_name}")
    try:
        return int(text)
    except ValueError as exc:
        raise ValueError(f"Invalid integer for {field_name}: {value}") from exc

def _build_sample_info(row: Dict[str, str], seen_ids: Set[str]) -> SampleInfo:
    sample_id = (row.get("sample_id") or "").strip()
    if not sample_id:
        raise ValueError("sample_id cannot be empty")
    if sample_id in seen_ids:
        raise ValueError(f"Duplicate sample_id detected: {sample_id}")
    
    run_mode = (row.get("run_mode") or "").strip().upper()
    if run_mode != "PE":
        raise ValueError(f"Unsupported run_mode '{run_mode}'. Currently only 'PE' is supported.")
    
    fastq1_raw = (row.get("path_FASTQ") or "").strip()
    if not fastq1_raw:
        raise ValueError("path_FASTQ cannot be empty")
    fastq1_path = Path(fastq1_raw).expanduser()
    if not fastq1_path.exists():
        raise ValueError(f"FASTQ file does not exist: {fastq1_path}")
    fastq1_path = fastq1_path.resolve()
    fastq2_path = infer_fastq2_path(fastq1_path).resolve()
    if not fastq2_path.exists():
        raise ValueError(f"Inferred FASTQ2 file does not exist: {fastq2_path}")
    
    group_value = (row.get("group") or "").strip() or "NA"
    tech_rep = _parse_int_field(row.get("tech_rep"), "tech_rep")
    bio_rep = _parse_int_field(row.get("bio_rep"), "bio_rep")
    batch_seq = _parse_int_field(row.get("batch_seq"), "batch_seq")
    
    seen_ids.add(sample_id)
    
    return SampleInfo(
        sample_id=sample_id,
        run_mode=run_mode,
        fastq1=fastq1_path,
        fastq2=fastq2_path,
        group=group_value,
        tech_rep=tech_rep,
        bio_rep=bio_rep,
        batch_seq=batch_seq
    )

def load_sample_info(sample_info_path: Path) -> Tuple[List[SampleInfo], List[str]]:
    """
    Load and validate sample information entries from a CSV file.
    Returns a tuple of (valid_samples, skipped_row_messages)
    """
    valid_samples: List[SampleInfo] = []
    skipped_rows: List[str] = []
    seen_ids: Set[str] = set()
    
    with sample_info_path.open("r", newline='') as csv_file:
        reader = csv.DictReader(csv_file)
        if reader.fieldnames is None:
            raise ValueError("Sample info CSV is missing a header row.")
        normalized_headers = [header.strip() if header else header for header in reader.fieldnames]
        reader.fieldnames = normalized_headers
        missing_columns = SAMPLE_INFO_REQUIRED_COLUMNS - set(normalized_headers)
        if missing_columns:
            missing = ", ".join(sorted(missing_columns))
            raise ValueError(f"Sample info CSV is missing required columns: {missing}")
        
        for row_index, row in enumerate(reader, start=2):
            try:
                sample_info = _build_sample_info(row, seen_ids)
                valid_samples.append(sample_info)
            except ValueError as exc:
                skipped_rows.append(f"Row {row_index}: {exc}")
    
    return valid_samples, skipped_rows

def create_sample_info_csv(
    sample_entries: List[SampleInfo],
    output_path: Path
) -> bool:
    """
    Create a Sample_info.csv file from sample entries for statistical analysis.
    
    Args:
        sample_entries: List of SampleInfo objects
        output_path: Path where to write the CSV file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline='') as csv_file:
            writer = csv.DictWriter(
                csv_file,
                fieldnames=["sample_id", "group", "tech_rep", "bio_rep", "batch_seq"]
            )
            writer.writeheader()
            for entry in sample_entries:
                writer.writerow({
                    "sample_id": entry.sample_id,
                    "group": entry.group,
                    "tech_rep": entry.tech_rep,
                    "bio_rep": entry.bio_rep,
                    "batch_seq": entry.batch_seq
                })
        logger.info(f"Created Sample_info.csv at {output_path}")
        return True
    except Exception as exc:
        logger.error(f"Failed to create Sample_info.csv at {output_path}: {exc}")
        return False

def merge_count_tables(
    sample_entries: List[SampleInfo],
    output_dir: Path,
    output_path: Optional[Path] = None
) -> Optional[Path]:
    """
    Merge individual per-sample count tables into a single wide-format table
    where rows are variants and columns include raw counts and sample-provided RPM values.
    """
    if not sample_entries:
        logger.warning("No sample entries were provided for count merging.")
        return None
    
    merged_records: Dict[Tuple[str, str], Dict[str, Any]] = {}
    sample_order: List[str] = []
    
    for entry in sample_entries:
        sample_order.append(entry.sample_id)
        counts_path = output_dir / entry.sample_id / f"{entry.sample_id}.counts.tsv"
        if not counts_path.exists():
            logger.warning(f"Count file not found for sample {entry.sample_id}: {counts_path}. Skipping.")
            continue
        try:
            with counts_path.open("r", newline='') as count_file:
                reader = csv.DictReader(count_file)
                required_cols = {"ID_WLG", "peptide", "count", "RPM"}
                if not reader.fieldnames or not required_cols.issubset(set(reader.fieldnames)):
                    logger.error(f"Count file for sample {entry.sample_id} is missing required columns. Skipping.")
                    continue
                for row in reader:
                    id_wlg = (row.get("ID_WLG") or "").strip()
                    if not id_wlg or id_wlg.lower() == "unassigned":
                        continue
                    key = (id_wlg, row["peptide"])
                    record = merged_records.setdefault(
                        key,
                        {
                            "ID_WLG": id_wlg,
                            "peptide": row["peptide"],
                            "counts": {},
                            "rpms": {}
                        }
                    )
                    try:
                        record["counts"][entry.sample_id] = int(row.get("count", 0) or 0)
                    except ValueError:
                        logger.warning(
                            f"Invalid count value in {counts_path} for variant {row['ID_WLG']}. Treating as 0."
                        )
                        record["counts"][entry.sample_id] = 0
                    try:
                        record["rpms"][entry.sample_id] = float(row.get("RPM", 0) or 0)
                    except ValueError:
                        logger.warning(
                            f"Invalid RPM value in {counts_path} for variant {row['ID_WLG']}. Treating as 0."
                        )
                        record["rpms"][entry.sample_id] = 0.0
        except Exception as exc:
            logger.error(f"Failed to read count file for sample {entry.sample_id}: {exc}")
            continue
    
    if not merged_records:
        logger.error("No count files could be read; merged counts table will not be created.")
        return None
    
    output_path = output_path or (output_dir / "merged.counts.tsv")
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline='') as merged_file:
            fieldnames = ["ID_WLG", "peptide"]
            rpm_columns: List[str] = []
            for sample_id in sample_order:
                fieldnames.append(sample_id)
                rpm_columns.append(f"{sample_id}_RPM")
            fieldnames.extend(rpm_columns)
            writer = csv.DictWriter(merged_file, fieldnames=fieldnames)
            writer.writeheader()
            for key in sorted(merged_records.keys()):
                record = merged_records[key]
                row = {
                    "ID_WLG": record["ID_WLG"],
                    "peptide": record["peptide"]
                }
                for sample_id in sample_order:
                    row[sample_id] = record["counts"].get(sample_id, 0)
                    row[f"{sample_id}_RPM"] = record["rpms"].get(sample_id, 0.0)
                writer.writerow(row)
        logger.info(f"Merged counts table written to {output_path}")
        return output_path
    except Exception as exc:
        logger.error(f"Failed to write merged counts table to {output_path}: {exc}")
        return None

def setup_logging(
    logger_instance: logging.Logger,
    level_str: str = "INFO",
    log_file: Optional[Path] = None
) -> Optional[logging.Handler]:
    """
    Configures the logging for the FastQ processing pipeline.
    
    Args:
        logger_instance (logging.Logger): The logger instance to configure.
        level_str (str): The logging level as a string (e.g., "DEBUG", "INFO").
        log_file (Path or None): Optional path to a file where logs should be written.
    """
    # Set the logger level based on the provided string
    log_level = getattr(logging, level_str.upper(), logging.INFO)
    logger_instance.setLevel(log_level)

    formatter = logging.Formatter(LOG_FORMAT)

    # Console Handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger_instance.addHandler(console_handler)

    file_handler: Optional[logging.Handler] = None

    # File Handler
    if log_file:
        try:
            log_file.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(log_file, mode='w')
            file_handler.setFormatter(formatter)
            logger_instance.addHandler(file_handler)
            logger_instance.info(f"Logging to console and to file: {log_file}")
        except Exception as e:
            logger_instance.error(f"Failed to set up file logging to {log_file}: {e}", exc_info=False)
            logger_instance.info("Continuing with console logging only.")
    else:
        logger_instance.info("Logging to console only.")

    return file_handler

def load_config(config_path=None):
    """Load configuration from file or use default."""
    config = None  # Initialize config variable
    
    if config_path and os.path.exists(config_path):
        try:
            with open(config_path, 'r') as f:
                config = json.load(f)
            logger.debug(f"Loaded configuration from {config_path}")
        except Exception as exc:
            logger.error(f"Error loading config file {config_path}: {exc}", exc_info=True)
            sys.exit(1)
    else:
        # No config path provided or file does not exist; use default config from package data
        try:
            with pkg_resources.open_text("capscreen", "config.json") as f:
                config = json.load(f)
            logger.debug("Loaded default configuration from package data.")
            
            # Try to use the config_path from the loaded config
            if config and "config_path" in config and os.path.exists(config["config_path"]):
                try:
                    with open(config["config_path"], 'r') as f:
                        config = json.load(f)
                    logger.debug(f"Loaded configuration from {config['config_path']}")
                except Exception as exc:
                    logger.warning(f"Could not load config from {config['config_path']}, using package default: {exc}")
        except Exception as exc:
            logger.error("Error: Default config file not found in package data.", exc_info=True)
            sys.exit(1)
    
    if config is None:
        logger.error("Failed to load any configuration.")
        sys.exit(1)
        
    return config
def sample_results_exist(output_dir: Path, sample_name: str, config: Optional[Dict[str, Any]] = None) -> bool:
    """
    Check whether the final count table for a sample already exists.
    Also checks for SAM file and PEAR merge output to avoid rerunning completed steps.
    """
    if not sample_name:
        return False
    sample_dir = output_dir / sample_name
    
    # Check for final count table - if exists, skip entire pipeline
    counts_file = sample_dir / f"{sample_name}.counts.tsv"
    if counts_file.exists() and counts_file.stat().st_size > 0:
        logger.info(
            "Detected existing results for sample %s at %s. Skipping re-run.",
            sample_name,
            counts_file
        )
        return True
    
    # Check for SAM file - if exists, skip alignment step (but will still run counting)
    sam_file = sample_dir / f"{sample_name}.sorted.sam"
    if sam_file.exists() and sam_file.stat().st_size > 0:
        logger.info(
            "Detected existing SAM file for sample %s at %s. Will skip alignment step.",
            sample_name,
            sam_file
        )
        # Don't return True here - we still want to run counting if count.tsv doesn't exist
    
    # Check for PEAR merge output - if exists, skip PEAR step (but will still run alignment)
    if config:
        pear_config = config.get('pear', {})
        pear_prefix = pear_config.get('output_prefix', 'pear')
        pear_output = sample_dir / f"{pear_prefix}.assembled.fastq.gz"
        if pear_output.exists() and pear_output.stat().st_size > 0:
            logger.info(
                "Detected existing PEAR merge output for sample %s at %s. Will skip PEAR step.",
                sample_name,
                pear_output
            )
            # Don't return True here - we still want to run alignment if SAM doesn't exist
    
    return False

def execute_full_pipeline(
    sample_name: str,
    fastq1: Path,
    fastq2: Path,
    output_dir: Path,
    reference_file: Path,
    config: Dict[str, Any],
    threads: Optional[int],
    keep_intermediate: bool,
    log_file: Optional[Path],
    cut_umi: bool = False
) -> bool:
    """
    Run the complete pipeline (QC, alignment, counting, reporting) for a single sample.
    Skips alignment step if SAM file already exists.
    """
    sample_dir = output_dir / sample_name
    sorted_sam = sample_dir / f"{sample_name}.sorted.sam"
    
    # Check if SAM file already exists - if so, skip alignment step
    if sorted_sam.exists() and sorted_sam.stat().st_size > 0:
        logger.info(
            "Detected existing SAM file for sample %s at %s. Skipping alignment step.",
            sample_name,
            sorted_sam
        )
    else:
        # Run QC and alignment steps
        sorted_sam = alignment_module.run_qc_and_alignment(
            fastq1,
            fastq2,
            output_dir,
            sample_name,
            config,
            threads,
            cut_umi=cut_umi,
        )
        if not sorted_sam:
            logger.error(f"QC/Alignment failed for sample {sample_name}.")
            return False
    
    output_file = sample_dir / f"{sample_name}.counts.tsv"
    
    # Use threads from parameter if provided, otherwise use config value
    # count_module.main will handle config fallback if n_threads is None
    n_threads = threads if threads is not None else config.get('threads')
    
    try:
        count_module.main(
            sorted_sam,
            reference_file,
            config,
            output_file,
            log_file=log_file,
            logger=logger,
            n_threads=n_threads
        )
    except Exception as e:
        logger.error(f"Counting failed for sample {sample_name}: {e}")
        return False
    
    try:
        logger.info("Generating HTML report...")
        generate_report_module.generate_report(str(sample_dir))
        logger.info("HTML report generated.")
    except Exception as e:
        logger.error(f"Report generation failed for sample {sample_name}: {e}", exc_info=True)
    
    if not keep_intermediate:
        alignment_module.cleanup_intermediate_files(sample_dir, sample_name)
        
    return True

def main():
    parser = argparse.ArgumentParser(
        description="CapScreen: QC, alignment, and variant counting pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--version", action="version", version=f"CapScreen v{__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Subcommand to run")

    # Parent parser for shared arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("--output-dir", type=Path, required=True, help="Base output directory")
    parent_parser.add_argument("--sample-name", required=False, help="Name of the sample")
    parent_parser.add_argument("--config", type=Path, default=Path("config.json"), help="Path to configuration file")
    parent_parser.add_argument("--threads", type=int, help="Number of threads to use")
    parent_parser.add_argument("--log-level", default="INFO", 
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level")

    # pipeline: full pipeline (QC, align, count)
    pipeline_parser = subparsers.add_parser("pipeline", parents=[parent_parser], help="Run full pipeline: QC, alignment, and count")
    pipeline_parser.add_argument("--fastq1", type=Path, required=False, help="Path to first FASTQ file")
    pipeline_parser.add_argument("--fastq2", type=Path, required=False, help="Path to second FASTQ file")
    pipeline_parser.add_argument("--sample-info", type=Path, required=False, help="CSV file describing multiple samples to process")
    pipeline_parser.add_argument("--merged-counts-output", type=Path, help="Path to write merged counts table when using --sample-info")
    pipeline_parser.add_argument("--reference-file", type=Path, required=True, help="Path to reference library file (CSV)")
    pipeline_parser.add_argument("--keep_intermediate", action="store_true", help="Keep intermediate FASTQ files (default: remove)")
    pipeline_parser.add_argument(
        "--cut-umi",
        action="store_true",
        help=(
            "After PEAR merging, run Cutadapt to trim reads to the variable region "
            "between configured flanking sequences and re-attach the flanks before alignment."
        ),
    )
    pipeline_parser.add_argument(
        "--stat",
        action="store_true",
        help="Run statistical analysis after creating merged counts table (requires --sample-info)"
    )

    # align: only QC and alignment
    align_parser = subparsers.add_parser("align", parents=[parent_parser], help="Run QC and alignment only")
    align_parser.add_argument("--fastq1", type=Path, required=True, help="Path to first FASTQ file")
    align_parser.add_argument("--fastq2", type=Path, required=True, help="Path to second FASTQ file")
    align_parser.add_argument(
        "--cut-umi",
        action="store_true",
        help=(
            "After PEAR merging, run Cutadapt to trim reads to the variable region "
            "between configured flanking sequences and re-attach the flanks before alignment."
        ),
    )

    # count: only counting
    count_parser = subparsers.add_parser("count", parents=[parent_parser], help="Run variant counting only")
    count_parser.add_argument("--sam-file", type=Path, required=True, help="Path to sorted SAM file")
    count_parser.add_argument("--reference-file", type=Path, required=True, help="Path to reference library file (CSV)")
    count_parser.add_argument("--output", type=Path, help="Path to output file (TSV)")

    # report: HTML report only
    report_parser = subparsers.add_parser("report", parents=[parent_parser], help="Generate HTML report from an existing sample directory")
    report_parser.add_argument("--output", type=Path, help="Output HTML file name (defaults to <sample>_report.html)")

    # stat: statistical analysis
    stat_parser = subparsers.add_parser("stat", parents=[parent_parser], help="Run statistical analysis on merged counts")
    stat_parser.add_argument("--counts", type=Path, required=True, help="Path to merged.counts.tsv file")
    stat_parser.add_argument("--sample-info", type=Path, required=True, help="Path to Sample_info.csv file")
    stat_parser.add_argument("--n-cpus", type=int, help="Number of CPUs to use for DESeq2")

    def normalize_path(path_value: Optional[Path]) -> Optional[Path]:
        if path_value is None:
            return None
        return path_value.expanduser().resolve()

    args = parser.parse_args()

    args.output_dir = normalize_path(args.output_dir)
    if isinstance(args.config, Path):
        args.config = normalize_path(args.config)

    for attr in ("fastq1", "fastq2", "reference_file", "sam_file", "output", "sample_info", "merged_counts_output", "counts"):
        if hasattr(args, attr):
            value = getattr(args, attr)
            if isinstance(value, Path):
                setattr(args, attr, normalize_path(value))

    # Command-specific validation for sample_name
    parser_map = {
        "align": align_parser,
        "count": count_parser,
        "report": report_parser
    }
    if args.command in parser_map and not args.sample_name:
        parser_map[args.command].error("--sample-name is required for this command.")

    if args.command == "pipeline":
        if args.sample_info:
            if not args.sample_info.exists():
                pipeline_parser.error(f"Sample info CSV not found: {args.sample_info}")
        else:
            missing_flags = []
            if not args.sample_name:
                missing_flags.append("--sample-name")
            if not args.fastq1:
                missing_flags.append("--fastq1")
            if not args.fastq2:
                missing_flags.append("--fastq2")
            if missing_flags:
                missing_str = ", ".join(missing_flags)
                pipeline_parser.error(f"{missing_str} must be provided when --sample-info is not used.")

    json_config = load_config(args.config)

    # Set up unified log file path
    batch_log_path: Optional[Path] = None
    if args.command == "pipeline" and args.sample_info:
        batch_label = args.sample_name if args.sample_name else "sample_info_batch"
        log_file = args.output_dir / f"{batch_label}.pipeline.log"
        batch_log_path = log_file
    else:
        log_identifier = args.sample_name if args.sample_name else "capscreen"
        log_file = args.output_dir / log_identifier / f"{log_identifier}.pipeline.log"

    log_file_handler = setup_logging(logger, args.log_level, log_file=log_file)

    if args.command == "pipeline":
        keep_intermediate = determine_keep_intermediate(getattr(args, 'keep_intermediate', False), json_config)
        if args.sample_info:
            logger.info(f"Sample info CSV provided: {args.sample_info}")
            try:
                sample_entries, skipped_rows = load_sample_info(args.sample_info)
            except ValueError as exc:
                logger.error(f"Failed to parse sample info CSV: {exc}")
                sys.exit(1)
            
            if skipped_rows:
                logger.warning(f"Skipping {len(skipped_rows)} row(s) from sample info CSV due to validation issues:")
                for row_msg in skipped_rows:
                    logger.warning(f"  - {row_msg}")
            
            if not sample_entries:
                logger.error("No valid rows were found in the provided sample info CSV. Aborting.")
                sys.exit(1)
            
            logger.info(f"Validated {len(sample_entries)} sample(s) for batch processing.")
            for entry in sample_entries:
                logger.info(
                    "Sample accepted: sample_id=%s run_mode=%s group=%s tech_rep=%d bio_rep=%d batch_seq=%d fastq1=%s fastq2=%s",
                    entry.sample_id,
                    entry.run_mode,
                    entry.group,
                    entry.tech_rep,
                    entry.bio_rep,
                    entry.batch_seq,
                    entry.fastq1,
                    entry.fastq2
                )

            if batch_log_path:
                detach_log_file_handler(logger, log_file_handler)
                log_file_handler = None
            
            skipped_existing = 0
            completed_samples = 0
            for entry in sample_entries:
                sample_dir = args.output_dir / entry.sample_id
                sample_dir.mkdir(parents=True, exist_ok=True)
                sample_log_file = sample_dir / f"{entry.sample_id}.pipeline.log"
                logger.info(f"Starting pipeline for sample {entry.sample_id}")
                if sample_results_exist(args.output_dir, entry.sample_id, json_config):
                    skipped_existing += 1
                    append_sample_log_to_batch(
                        batch_log_path,
                        sample_log_file,
                        entry.sample_id,
                        fallback_message="Sample skipped because existing results were detected.",
                        copy_sample_log=False
                    )
                    continue
                with sample_file_log_handler(sample_log_file):
                    success = execute_full_pipeline(
                        entry.sample_id,
                        entry.fastq1,
                        entry.fastq2,
                        args.output_dir,
                        args.reference_file,
                        json_config,
                        args.threads,
                        keep_intermediate,
                        sample_log_file,
                        cut_umi=getattr(args, "cut_umi", False),
                    )
                append_sample_log_to_batch(batch_log_path, sample_log_file, entry.sample_id)
                if not success:
                    logger.warning(f"Pipeline reported a failure for sample {entry.sample_id}. Checking for count file before proceeding.")
                count_file = args.output_dir / entry.sample_id / f"{entry.sample_id}.counts.tsv"
                if not count_file.exists():
                    logger.error(f"Count file not found for sample {entry.sample_id} at {count_file}. Moving on to next sample.")
                    continue
                completed_samples += 1
            logger.info(
                "Pipeline completed successfully for %d sample(s); skipped %d sample(s) with existing results.",
                completed_samples,
                skipped_existing
            )
            merged_output_path = args.merged_counts_output
            if merged_output_path is None:
                merged_output_path = args.output_dir / "merged.counts.tsv"
            merged_path = merge_count_tables(sample_entries, args.output_dir, merged_output_path)
            if merged_path:
                logger.info(f"Merged counts table saved to {merged_path}")
            else:
                logger.warning("Merged counts table could not be created.")
            
            # Track statistical analysis success status
            stat_success = None  # None = not attempted, True = succeeded, False = failed
            
            # Run statistical analysis if requested
            if getattr(args, 'stat', False):
                if not merged_path:
                    logger.error("Cannot run statistical analysis: merged counts table was not created.")
                    stat_success = False
                else:
                    logger.info("Starting QC analysis before statistical analysis...")
                    # Create output directory for stat results (QC plots will be saved here too)
                    stat_output_dir = args.output_dir / "statistical_analysis"
                    stat_output_dir.mkdir(parents=True, exist_ok=True)
                    
                    # Run QC analysis on reports in output_dir, save plots to stat_output_dir
                    qc_success = qc_module.run_qc_analysis(args.output_dir, stat_output_dir, logger_instance=logger)
                    if qc_success:
                        logger.info("QC analysis completed successfully")
                    else:
                        logger.warning("QC analysis completed with warnings or errors, continuing with statistical analysis...")
                    
                    logger.info("Starting statistical analysis...")
                    # Use the original Sample_info.csv file provided via --sample-info (preserve all columns)
                    # Do NOT create a new file that would overwrite the original
                    if args.sample_info and args.sample_info.exists():
                        sample_info_path = args.sample_info
                        logger.info(f"Using original Sample_info.csv at {sample_info_path} (all columns preserved)")
                        
                        # Run statistical analysis using stat environment
                        try:
                            success = run_statistical_analysis_with_stat_env(
                                counts_file=merged_path,
                                sample_info_file=sample_info_path,
                                output_dir=stat_output_dir,
                                config_file=args.config,
                                n_cpus=args.threads,
                                verbose=(args.log_level == "DEBUG"),
                                logger_instance=logger
                            )
                            stat_success = success
                            if success:
                                logger.info(f"Statistical analysis completed. Results saved to {stat_output_dir}")
                            else:
                                logger.warning("Statistical analysis completed with warnings or errors.")
                        except Exception as e:
                            logger.error(f"Statistical analysis failed: {e}", exc_info=True)
                            stat_success = False
                    else:
                        logger.error("Original Sample_info.csv file not found. Cannot proceed with statistical analysis.")
                        stat_success = False
            
            # Determine statistical analysis status for summary
            if stat_success is True:
                stat_status = "Completed"
            elif stat_success is False:
                stat_status = "Failed"
            elif getattr(args, 'stat', False):
                stat_status = "Skipped (merged counts table not available)"
            else:
                stat_status = "Skipped"
            
            append_batch_summary(
                batch_log_path,
                [
                    f"Completed samples: {completed_samples}",
                    f"Skipped samples (existing results): {skipped_existing}",
                    f"Merged counts table: {merged_path}" if merged_path else "Merged counts table was not created.",
                    f"Statistical analysis: {stat_status}"
                ]
            )
            sys.exit(0)
        else:
            if sample_results_exist(args.output_dir, args.sample_name, json_config):
                logger.info("To re-run the pipeline for this sample, remove the existing results or move them elsewhere.")
                sys.exit(0)
            success = execute_full_pipeline(
                args.sample_name,
                args.fastq1,
                args.fastq2,
                args.output_dir,
                args.reference_file,
                json_config,
                args.threads,
                keep_intermediate,
                log_file,
                cut_umi=getattr(args, "cut_umi", False),
            )
            sys.exit(0 if success else 1)
    elif args.command == "align":
        sorted_sam = alignment_module.run_qc_and_alignment(
            args.fastq1,
            args.fastq2,
            args.output_dir,
            args.sample_name,
            json_config,
            args.threads,
            cut_umi=getattr(args, "cut_umi", False),
        )
        if not sorted_sam:
            logger.error("QC/Alignment failed.")
            sys.exit(1)
        logger.info(f"QC and alignment completed successfully. Sorted SAM: {sorted_sam}")
        sys.exit(0)
    elif args.command == "count":
        # Use provided output file or default
        output_file = args.output if args.output else args.sam_file.parent / f"{args.sam_file.stem}.counts.tsv"
        # Use --threads if provided, otherwise use config value (count_module.main will handle config fallback)
        n_threads = args.threads if args.threads is not None else json_config.get('threads')
        try:
            count_module.main(
                args.sam_file,
                args.reference_file,
                json_config,
                output_file,
                log_file=log_file,
                logger=logger,
                n_threads=n_threads
            )
        except Exception as e:
            logger.error(f"Counting failed: {e}")
            sys.exit(1)
        logger.info("Counting completed successfully.")
        sys.exit(0)
    elif args.command == "report":
        sample_dir = args.output_dir / args.sample_name
        if not sample_dir.exists():
            logger.error(f"Sample directory not found: {sample_dir}")
            sys.exit(1)
        try:
            logger.info("Generating HTML report...")
            output_arg = str(args.output) if args.output else None
            generate_report_module.generate_report(str(sample_dir), output=output_arg)
            logger.info("HTML report generated.")
        except Exception as e:
            logger.error(f"Report generation failed: {e}", exc_info=True)
            sys.exit(1)
        sys.exit(0)
    elif args.command == "stat":
        # Validate inputs
        if not args.counts.exists():
            stat_parser.error(f"Counts file not found: {args.counts}")
        if not args.sample_info.exists():
            stat_parser.error(f"Sample info file not found: {args.sample_info}")
        
        # Create output directory
        args.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run statistical analysis using stat environment
        try:
            n_cpus_val = getattr(args, 'n_cpus', None) or args.threads
            success = run_statistical_analysis_with_stat_env(
                counts_file=args.counts,
                sample_info_file=args.sample_info,
                output_dir=args.output_dir,
                config_file=args.config,
                n_cpus=n_cpus_val,
                verbose=(args.log_level == "DEBUG"),
                logger_instance=logger
            )
            if success:
                logger.info("Statistical analysis completed successfully.")
                sys.exit(0)
            else:
                logger.error("Statistical analysis failed.")
                sys.exit(1)
        except Exception as e:
            logger.error(f"Statistical analysis failed: {e}", exc_info=True)
            sys.exit(1)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
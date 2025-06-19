#!/usr/bin/env python3
import argparse
import json
import logging
import os
import sys
import subprocess
import shlex
from pathlib import Path
from typing import Optional, Dict, Any, Tuple
import importlib.resources as pkg_resources
from capscreen.version import __version__

# Global Logger
logger = logging.getLogger("FastQProcessor")

def setup_logging(logger_instance: logging.Logger, level_str: str = "INFO", log_file: Optional[Path] = None) -> None:
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

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] [%(name)s:%(lineno)d] %(message)s")

    # Console Handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger_instance.addHandler(console_handler)

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

def run_command(command: list, step_name: str, log_file: Optional[Path] = None) -> bool:
    """
    Run a command and log its output.
    
    Args:
        command (list): Command to run as a list of strings
        step_name (str): Name of the processing step
        log_file (Path, optional): Path to log file for command output
        
    Returns:
        bool: True if command succeeded, False otherwise
    """
    logger.info(f"[{step_name}] Running command: {' '.join(shlex.quote(str(s)) for s in command)}")
    
    try:
        process = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True
        )
        
        if process.stdout:
            logger.debug(f"[{step_name}] STDOUT:\n{process.stdout.strip()}")
        if process.stderr:
            logger.debug(f"[{step_name}] STDERR:\n{process.stderr.strip()}")
            
        logger.info(f"[{step_name}] Successfully completed.")
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"[{step_name}] Command failed with exit code {e.returncode}.")
        if e.stdout:
            logger.error(f"[{step_name}] STDOUT:\n{e.stdout.strip()}")
        if e.stderr:
            logger.error(f"[{step_name}] STDERR:\n{e.stderr.strip()}")
        return False
    except Exception as e:
        logger.error(f"[{step_name}] An unexpected error occurred: {e}", exc_info=True)
        return False

def run_fastp(fastq1: Path, fastq2: Path, sample_dir: Path, config: Dict[str, Any], threads: int) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Run FASTP on the input FASTQ files.
    
    Args:
        fastq1 (Path): Path to first FASTQ file
        fastq2 (Path): Path to second FASTQ file
        sample_dir (Path): Sample-specific output directory
        config (Dict): Configuration dictionary
        threads (int): Number of threads to use
        
    Returns:
        Tuple[Optional[Path], Optional[Path]]: Paths to trimmed FASTQ files
    """
    fastp_config = config['fastp']
    prefix = fastp_config['output_prefix']
    
    r1_trimmed = sample_dir / f"{prefix}_R1.fastq.gz"
    r2_trimmed = sample_dir / f"{prefix}_R2.fastq.gz"
    
    fastp_cmd = [
        "fastp",
        "--thread", str(threads),
        "--in1", str(fastq1),
        "--in2", str(fastq2),
        "--out1", str(r1_trimmed),
        "--out2", str(r2_trimmed),
        "--compression", str(fastp_config['compression']),
        "--length_required", str(fastp_config['length_required']),
        "--dup_calc_accuracy", str(fastp_config['dup_calc_accuracy'])
    ]
    
    if fastp_config['html_report']:
        fastp_cmd.extend(["--html", str(sample_dir / f"{prefix}.html")])
    if fastp_config['json_report']:
        fastp_cmd.extend(["--json", str(sample_dir / f"{prefix}.json")])
    
    if run_command(fastp_cmd, "FASTP"):
        return r1_trimmed, r2_trimmed
    return None, None

def run_pear(r1_trimmed: Path, r2_trimmed: Path, sample_dir: Path, config: Dict[str, Any], threads: int) -> Optional[Path]:
    """
    Run PEAR to merge the trimmed FASTQ files.

    Args:
        r1_trimmed (Path): Path to trimmed R1 FASTQ
        r2_trimmed (Path): Path to trimmed R2 FASTQ
        sample_dir (Path): Sample-specific output directory
        config (Dict): Configuration dictionary
        threads (int): Number of threads to use

    Returns:
        Optional[Path]: Path to assembled FASTQ file
    """
    pear_config = config['pear']
    prefix = pear_config['output_prefix']
    pear_base = sample_dir / prefix
    log_file = sample_dir / f"{prefix}.pear.log"
    
    # PEAR command
    pear_cmd = [
        "pear",
        "-f", str(r1_trimmed),
        "-r", str(r2_trimmed),
        "-o", str(pear_base),
        "-j", str(threads)
    ]
    
    try:
        logger.info("Running PEAR")
        with open(log_file, 'w') as log:
            process = subprocess.run(
                pear_cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True
            )
            
            if process.returncode != 0:
                logger.error(f"PEAR failed with return code {process.returncode}")
                return None
            
        # Compress all PEAR output files
        files_to_compress = [
            f"{pear_base}.assembled.fastq",
            f"{pear_base}.discarded.fastq",
            f"{pear_base}.unassembled.forward.fastq",
            f"{pear_base}.unassembled.reverse.fastq"
        ]
        
        for file in files_to_compress:
            if Path(file).exists():
                gzip_cmd = ["gzip", "-f", file]  # -f to force overwrite
                if not run_command(gzip_cmd, f"Gzip {file}"):
                    logger.warning(f"Failed to compress {file}")
        
        # Return path to compressed assembled file
        assembled_fastq_gz = f"{pear_base}.assembled.fastq.gz"
        if Path(assembled_fastq_gz).exists():
            return Path(assembled_fastq_gz)
        else:
            logger.error("Compressed assembled file not found")
            return None
            
    except Exception as e:
        logger.error(f"Error in PEAR pipeline: {e}", exc_info=True)
        return None

def validate_config(config: Dict[str, Any]) -> bool:
    """
    Validate the configuration dictionary structure.
    
    Args:
        config (Dict[str, Any]): Configuration dictionary to validate
        
    Returns:
        bool: True if configuration is valid, False otherwise
    """
    required_sections = ['fastp', 'pear', 'bowtie2', 'samtools', 'reference']
    required_fastp_keys = ['output_prefix', 'compression', 'length_required', 'dup_calc_accuracy']
    required_bowtie2_keys = ['output_prefix', 'mode', 'sensitivity', 'np', 'n_ceil']
    required_samtools_keys = ['sort_memory']
    required_reference_keys = ['genome']
    
    try:
        # Check required sections
        for section in required_sections:
            if section not in config:
                logger.error(f"Missing required configuration section: {section}")
                return False
        
        # Check required keys in each section
        for key in required_fastp_keys:
            if key not in config['fastp']:
                logger.error(f"Missing required key in fastp section: {key}")
                return False
        
        for key in required_bowtie2_keys:
            if key not in config['bowtie2']:
                logger.error(f"Missing required key in bowtie2 section: {key}")
                return False
        
        for key in required_samtools_keys:
            if key not in config['samtools']:
                logger.error(f"Missing required key in samtools section: {key}")
                return False
        
        for key in required_reference_keys:
            if key not in config['reference']:
                logger.error(f"Missing required key in reference section: {key}")
                return False
        
        return True
    except Exception as e:
        logger.error(f"Error validating configuration: {e}")
        return False

def validate_paths(fastq1: Path, fastq2: Path, output_dir: Path) -> bool:
    """
    Validate input and output paths.
    
    Args:
        fastq1 (Path): Path to first FASTQ file
        fastq2 (Path): Path to second FASTQ file
        output_dir (Path): Output directory path
        
    Returns:
        bool: True if all paths are valid, False otherwise
    """
    try:
        # Check input files
        if not fastq1.exists():
            logger.error(f"Input file does not exist: {fastq1}")
            return False
        if not fastq2.exists():
            logger.error(f"Input file does not exist: {fastq2}")
            return False
        
        # Check output directory
        if output_dir.exists() and not os.access(output_dir, os.W_OK):
            logger.error(f"No write permission for output directory: {output_dir}")
            return False
        
        return True
    except Exception as e:
        logger.error(f"Error validating paths: {e}")
        return False

def validate_threads(threads: int) -> int:
    """
    Validate and adjust thread count.
    
    Args:
        threads (int): Requested number of threads
        
    Returns:
        int: Validated thread count
    """
    try:
        threads = int(threads)
        if threads < 1:
            logger.warning("Thread count must be positive. Setting to 1.")
            return 1
        # Get system CPU count and limit threads
        cpu_count = os.cpu_count() or 1
        if threads > cpu_count:
            logger.warning(f"Thread count ({threads}) exceeds available CPUs ({cpu_count}). Setting to {cpu_count}.")
            return cpu_count
        return threads
    except (ValueError, TypeError):
        logger.warning("Invalid thread count. Setting to 1.")
        return 1

def run_bowtie2_and_samtools(assembled_fastq: Path, sample_dir: Path, config: Dict[str, Any], threads: int) -> bool:
    """
    Run Bowtie2 alignment and Samtools processing with piped commands.

    Args:
        assembled_fastq (Path): Path to assembled FASTQ file
        sample_dir (Path): Sample-specific output directory
        config (Dict): Configuration dictionary
        threads (int): Number of threads to use

    Returns:
        bool: True if processing succeeded, False otherwise
    """
    bowtie2_config = config['bowtie2']
    samtools_config = config['samtools']
    ref_config = config['reference']
    
    # Final output SAM file
    prefix = bowtie2_config['output_prefix']
    sorted_sam = sample_dir / f"{prefix}.sorted.sam"
    log_file = sample_dir / f"{prefix}.bowtie2.log"
    
    # Bowtie2 command with memory optimization
    bowtie2_cmd = [
        "bowtie2",
        "-x", ref_config['genome'],
        "-U", str(assembled_fastq),
        "--" + bowtie2_config['mode'],
        "--" + bowtie2_config['sensitivity'],
        "--np", str(bowtie2_config['np']),
        "--n-ceil", bowtie2_config['n_ceil'],
        "-p", str(threads)
    ]
    
    if bowtie2_config['xeq']:
        bowtie2_cmd.append("--xeq")
    if bowtie2_config['reorder']:
        bowtie2_cmd.append("--reorder")
    if bowtie2_config['N']:
        bowtie2_cmd.extend(["-N", str(bowtie2_config['N'])])
    if bowtie2_config['score_min']:
        bowtie2_cmd.extend(["--score-min", bowtie2_config['score_min']])
    
    # Samtools sort command for SAM
    samtools_sort_cmd = [
        "samtools", "sort",
        "-@", str(threads),
        "-m", samtools_config['sort_memory'],
        "-o", str(sorted_sam),
        "-"
    ]
    
    # Run the pipeline
    try:
        logger.info("Running Bowtie2 -> Samtools sort pipeline")
        
        with open(log_file, 'w') as log:
            with subprocess.Popen(bowtie2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as p1, \
                subprocess.Popen(samtools_sort_cmd, stdin=p1.stdout, stderr=subprocess.PIPE) as p2:
                
                # Get stderr from all processes
                bowtie2_stderr = p1.stderr.read().decode()
                sort_stderr = p2.stderr.read().decode()
                
                # Write Bowtie2 stderr to log file
                if bowtie2_stderr:
                    log.write("Bowtie2 output:\n")
                    log.write(bowtie2_stderr)
                    log.write("\n")
                
                # Write Samtools sort stderr to log file
                if sort_stderr:
                    log.write("Samtools sort output:\n")
                    log.write(sort_stderr)
                    log.write("\n")
                
                # Wait for all processes to complete
                p1.wait()
                p2.wait()
                
                # Check for errors
                if p1.returncode != 0:
                    logger.error(f"Bowtie2 failed with return code {p1.returncode}")
                    return False
                if p2.returncode != 0:
                    logger.error(f"Samtools sort failed with return code {p2.returncode}")
                    return False
        
        return True
        
    except Exception as e:
        logger.error(f"Error in Bowtie2/Samtools pipeline: {e}", exc_info=True)
        return False

def process_fastq_files(fastq1: Path, fastq2: Path, output_dir: Path, sample_name: str, config: Dict[str, Any], threads: int) -> bool:
    """
    Process the FASTQ files with the given parameters.
    
    Args:
        fastq1 (Path): Path to the first FASTQ file
        fastq2 (Path): Path to the second FASTQ file
        output_dir (Path): Base output directory
        sample_name (str): Name of the sample
        config (Dict): Configuration dictionary
        threads (int): Number of threads to use
        
    Returns:
        bool: True if processing was successful, False otherwise
    """
    try:
        # Validate configuration
        if not validate_config(config):
            return False
        
        # Validate paths
        if not validate_paths(fastq1, fastq2, output_dir):
            return False
        
        # Validate and adjust thread count
        threads = validate_threads(threads)
        
        # Create sample-specific directory
        sample_dir = output_dir / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up logging for this sample
        log_file = sample_dir / f"{sample_name}.log"
        setup_logging(logger, level_str=config.get('log_level', 'INFO'), log_file=log_file)
        
        logger.info(f"Processing sample: {sample_name}")
        logger.info(f"Output directory: {sample_dir}")
        logger.info(f"Using {threads} threads")
        
        # Step 1: Run FASTP
        logger.info("Step 1: Running FASTP")
        r1_trimmed, r2_trimmed = run_fastp(fastq1, fastq2, sample_dir, config, threads)
        if not r1_trimmed or not r2_trimmed:
            return False
        
        # Step 2: Run PEAR
        logger.info("Step 2: Running PEAR")
        if not run_pear(r1_trimmed, r2_trimmed, sample_dir, config, threads):
            return False
        
        # Step 3: Run Bowtie2 and Samtools
        logger.info("Step 3: Running Bowtie2 and Samtools")
        if not run_bowtie2_and_samtools(r1_trimmed, sample_dir, config, threads):
            return False
        
        logger.info("All processing steps completed successfully")
        return True
        
    except Exception as e:
        logger.error(f"Error processing FASTQ files: {e}", exc_info=True)
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Process FASTQ files with FASTP, PEAR, and Bowtie2",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Version argument
    parser.add_argument("--version", action="version",
                        version=f"CapScreen v{__version__}")
    
    # Required arguments
    parser.add_argument("--fastq1", type=Path, required=True, help="Path to first FASTQ file")
    parser.add_argument("--fastq2", type=Path, required=True, help="Path to second FASTQ file")
    parser.add_argument("--output-dir", type=Path, required=True, help="Base output directory")
    parser.add_argument("--sample-name", required=True, help="Name of the sample")
    
    # Optional arguments
    parser.add_argument("--threads", type=int, help="Number of threads to use")
    parser.add_argument("--config", type=Path, default=Path("config.json"), help="Path to configuration file")
    parser.add_argument("--log-level", default="INFO", 
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Set the logging level")
    
    # Parse all arguments
    args = parser.parse_args()
    
    # Load configuration from JSON
    json_config = load_config(args.config)
    
    # Set up logging
    setup_logging(logger, args.log_level)
    
    # Process FASTQ files
    if process_fastq_files(args.fastq1, args.fastq2, args.output_dir, args.sample_name, json_config, args.threads):
        logger.info("FASTQ processing completed successfully")
        sys.exit(0)
    else:
        logger.error("FASTQ processing failed")
        sys.exit(1)

if __name__ == "__main__":
    main()
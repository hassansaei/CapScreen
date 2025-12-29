#!/usr/bin/env python3
"""
Alignment and QC utilities used by the CapScreen CLI.
"""
import logging
import os
import shlex
import subprocess
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

logger = logging.getLogger("FastQProcessor")


def run_command(command: list, step_name: str, log_file: Optional[Path] = None) -> bool:
    """
    Run a shell command and log its output.
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

    except subprocess.CalledProcessError as err:
        logger.error(f"[{step_name}] Command failed with exit code {err.returncode}.")
        if err.stdout:
            logger.error(f"[{step_name}] STDOUT:\n{err.stdout.strip()}")
        if err.stderr:
            logger.error(f"[{step_name}] STDERR:\n{err.stderr.strip()}")
        return False
    except Exception as exc:
        logger.error(f"[{step_name}] An unexpected error occurred: {exc}", exc_info=True)
        return False


def run_fastp(
    fastq1: Path,
    fastq2: Path,
    sample_dir: Path,
    config: Dict[str, Any],
    threads: int
) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Run FASTP on the input FASTQ files.
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


def run_pear(
    r1_trimmed: Path,
    r2_trimmed: Path,
    sample_dir: Path,
    config: Dict[str, Any],
    threads: int
) -> Optional[Path]:
    """
    Run PEAR to merge the trimmed FASTQ files.
    """
    pear_config = config['pear']
    prefix = pear_config['output_prefix']
    pear_base = sample_dir / prefix
    log_file = sample_dir / f"{prefix}.pear.log"

    pear_cmd = [
        "pear",
        "-f", str(r1_trimmed),
        "-r", str(r2_trimmed),
        "-o", str(pear_base),
        "-j", str(threads)
    ]

    try:
        logger.info("Running PEAR")
        logger.info(f"[PEAR] Running command: {' '.join(shlex.quote(str(s)) for s in pear_cmd)}")
        with open(log_file, 'w') as log_handle:
            process = subprocess.run(
                pear_cmd,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True
            )

            if process.returncode != 0:
                logger.error(f"PEAR failed with return code {process.returncode}")
                return None

        files_to_compress = [
            f"{pear_base}.assembled.fastq",
            f"{pear_base}.discarded.fastq",
            f"{pear_base}.unassembled.forward.fastq",
            f"{pear_base}.unassembled.reverse.fastq"
        ]

        for file_path in files_to_compress:
            path_obj = Path(file_path)
            if path_obj.exists():
                gzip_cmd = ["gzip", "-f", file_path]
                if not run_command(gzip_cmd, f"Gzip {file_path}"):
                    logger.warning(f"Failed to compress {file_path}")

        assembled_fastq_gz = f"{pear_base}.assembled.fastq.gz"
        if Path(assembled_fastq_gz).exists():
            return Path(assembled_fastq_gz)

        logger.error("Compressed assembled file not found")
        return None

    except Exception as exc:
        logger.error(f"Error in PEAR pipeline: {exc}", exc_info=True)
        return None


def run_cutadapt_umi(
    assembled_fastq: Path,
    sample_dir: Path,
    config: Dict[str, Any],
    threads: int
) -> Optional[Path]:
    """
    Run Cutadapt to extract the variable region between flanking sequences
    and then re-attach the flanks to each read.

    Workflow:
      1. Use Cutadapt to trim everything outside the configured flanks
         on the PEAR-assembled reads (and discard reads without both flanks).
      2. Re-add the known flanks (from the config) to each trimmed read so
         downstream steps that expect flanking sequences still work.
    """
    flanking_cfg = config.get("flanking_sequences") or {}
    flank_5p = flanking_cfg.get("flank_5p")
    flank_3p = flanking_cfg.get("flank_3p")

    if not flank_5p or not flank_3p:
        logger.error(
            "Flanking sequences (flank_5p and flank_3p) must be defined in the "
            "configuration in order to use UMI cutting."
        )
        return None

    trimmed_fastq = sample_dir / "pear.trimmed.fastq"
    cutadapt_cmd = [
        "cutadapt",
        "-g", flank_5p,
        "-a", flank_3p,
        "--discard-untrimmed",
        "-o", str(trimmed_fastq),
        str(assembled_fastq),
    ]

    # Use multiple cores if available (Cutadapt uses -j/--cores)
    if threads and threads > 1:
        cutadapt_cmd[1:1] = ["-j", str(threads)]

    if not run_command(cutadapt_cmd, "CUTADAPT_UMI"):
        return None

    # Re-add flanks to each read so that downstream steps (alignment + counting)
    # still see the flanking regions around the variable region.
    restored_fastq = sample_dir / "pear.trimmed_with_flanks.fastq"
    try:
        logger.info("Re-attaching flanking sequences to Cutadapt-trimmed reads")
        with trimmed_fastq.open("r") as in_f, restored_fastq.open("w") as out_f:
            while True:
                header = in_f.readline()
                if not header:
                    break
                seq = in_f.readline()
                plus = in_f.readline()
                qual = in_f.readline()

                if not (header and seq and plus and qual):
                    logger.warning(
                        "Encountered incomplete FASTQ record while restoring flanks; "
                        "stopping early."
                    )
                    break

                seq = seq.rstrip("\n\r")
                qual = qual.rstrip("\n\r")

                new_seq = f"{flank_5p}{seq}{flank_3p}"
                flank_5p_qual = "I" * len(flank_5p)
                flank_3p_qual = "I" * len(flank_3p)
                new_qual = f"{flank_5p_qual}{qual}{flank_3p_qual}"

                out_f.write(header)
                out_f.write(new_seq + "\n")
                out_f.write(plus)
                out_f.write(new_qual + "\n")

    except Exception as exc:
        logger.error(f"Failed to restore flanks on Cutadapt-trimmed reads: {exc}", exc_info=True)
        return None

    # Compress the trimmed FASTQ file to save space
    if trimmed_fastq.exists():
        gzip_cmd = ["gzip", "-f", str(trimmed_fastq)]
        if not run_command(gzip_cmd, f"Gzip {trimmed_fastq}"):
            logger.warning(f"Failed to compress {trimmed_fastq}")

    # Compress the final restored FASTQ to save space and match other steps
    restored_fastq_gz = f"{restored_fastq}.gz"
    gzip_cmd = ["gzip", "-f", str(restored_fastq)]
    if not run_command(gzip_cmd, f"Gzip {restored_fastq}"):
        logger.warning(f"Failed to compress {restored_fastq}; using uncompressed file instead.")
        return restored_fastq

    return Path(restored_fastq_gz)


def validate_config(config: Dict[str, Any]) -> bool:
    """
    Validate the configuration dictionary structure.
    """
    required_sections = ['fastp', 'pear', 'bowtie2', 'samtools', 'reference']
    required_fastp_keys = ['output_prefix', 'compression', 'length_required', 'dup_calc_accuracy']
    required_bowtie2_keys = ['output_prefix', 'mode', 'sensitivity', 'np', 'n_ceil']
    required_samtools_keys = ['sort_memory']
    required_reference_keys = ['genome']

    try:
        for section in required_sections:
            if section not in config:
                logger.error(f"Missing required configuration section: {section}")
                return False

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
    except Exception as exc:
        logger.error(f"Error validating configuration: {exc}", exc_info=True)
        return False


def validate_paths(fastq1: Path, fastq2: Path, output_dir: Path) -> bool:
    """
    Validate input and output paths for alignment.
    """
    try:
        if not fastq1.exists():
            logger.error(f"Input file does not exist: {fastq1}")
            return False
        if not fastq2.exists():
            logger.error(f"Input file does not exist: {fastq2}")
            return False

        if output_dir.exists() and not os.access(output_dir, os.W_OK):
            logger.error(f"No write permission for output directory: {output_dir}")
            return False

        return True
    except Exception as exc:
        logger.error(f"Error validating paths: {exc}", exc_info=True)
        return False


def validate_threads(threads: Optional[int]) -> int:
    """
    Validate and adjust thread count.
    """
    try:
        threads = int(threads) if threads is not None else 1
        if threads < 1:
            logger.warning("Thread count must be positive. Setting to 1.")
            return 1
        cpu_count = os.cpu_count() or 1
        if threads > cpu_count:
            logger.warning(f"Thread count ({threads}) exceeds available CPUs ({cpu_count}). Setting to {cpu_count}.")
            return cpu_count
        return threads
    except (ValueError, TypeError):
        logger.warning("Invalid thread count. Setting to 1.")
        return 1


def run_bowtie2_and_samtools(
    assembled_fastq: Path,
    sample_dir: Path,
    sample_name: str,
    config: Dict[str, Any],
    threads: int
) -> bool:
    """
    Run Bowtie2 alignment and Samtools processing with piped commands.
    """
    bowtie2_config = config['bowtie2']
    samtools_config = config['samtools']
    ref_config = config['reference']

    sorted_sam = sample_dir / f"{sample_name}.sorted.sam"
    log_file = sample_dir / f"{sample_name}.bowtie2.log"

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

    samtools_sort_cmd = [
        "samtools", "sort",
        "-@", str(threads),
        "-m", samtools_config['sort_memory'],
        "-o", str(sorted_sam),
        "-"
    ]

    try:
        logger.info("Running Bowtie2 -> Samtools sort pipeline")
        bowtie2_str = ' '.join(shlex.quote(str(s)) for s in bowtie2_cmd)
        samtools_str = ' '.join(shlex.quote(str(s)) for s in samtools_sort_cmd)
        logger.info(f"[Bowtie2/Samtools] Running command: {bowtie2_str} | {samtools_str}")
        with open(log_file, 'w') as log_handle:
            with subprocess.Popen(bowtie2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as p1, \
                    subprocess.Popen(samtools_sort_cmd, stdin=p1.stdout, stderr=subprocess.PIPE) as p2:

                bowtie2_stderr = p1.stderr.read().decode()
                sort_stderr = p2.stderr.read().decode()

                if bowtie2_stderr:
                    log_handle.write("Bowtie2 output:\n")
                    log_handle.write(bowtie2_stderr)
                    log_handle.write("\n")

                if sort_stderr:
                    log_handle.write("Samtools sort output:\n")
                    log_handle.write(sort_stderr)
                    log_handle.write("\n")

                p1.wait()
                p2.wait()

                if p1.returncode != 0:
                    logger.error(f"Bowtie2 failed with return code {p1.returncode}")
                    return False
                if p2.returncode != 0:
                    logger.error(f"Samtools sort failed with return code {p2.returncode}")
                    return False

        return True

    except Exception as exc:
        logger.error(f"Error in Bowtie2/Samtools pipeline: {exc}", exc_info=True)
        return False


def run_qc_and_alignment(
    fastq1: Path,
    fastq2: Path,
    output_dir: Path,
    sample_name: str,
    config: Dict[str, Any],
    threads: Optional[int],
    cut_umi: bool = False
) -> Optional[Path]:
    """
    Run QC and alignment steps (FASTP, PEAR, Bowtie2/Samtools).
    Skips PEAR step if PEAR output already exists.
    Returns the path to the sorted SAM file if successful.
    """
    if not validate_config(config):
        return None
    if not validate_paths(fastq1, fastq2, output_dir):
        return None
    threads = validate_threads(threads)
    sample_dir = output_dir / sample_name
    sample_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Processing sample: {sample_name}")
    logger.info(f"Output directory: {sample_dir}")
    logger.info(f"Using {threads} threads")

    # Check if PEAR output already exists
    pear_config = config.get('pear', {})
    pear_prefix = pear_config.get('output_prefix', 'pear')
    pear_base = sample_dir / pear_prefix
    assembled_fastq = Path(f"{pear_base}.assembled.fastq.gz")
    
    if assembled_fastq.exists() and assembled_fastq.stat().st_size > 0:
        logger.info(
            "Detected existing PEAR merge output at %s. Skipping FASTP and PEAR steps.",
            assembled_fastq
        )
    else:
        logger.info("Step 1: Running FASTP")
        r1_trimmed, r2_trimmed = run_fastp(fastq1, fastq2, sample_dir, config, threads)
        if not r1_trimmed or not r2_trimmed:
            return None

        logger.info("Step 2: Running PEAR")
        assembled_fastq = run_pear(r1_trimmed, r2_trimmed, sample_dir, config, threads)
        if not assembled_fastq:
            return None

    if cut_umi:
        logger.info("Step 3: Running Cutadapt-based UMI trimming")
        umi_fastq = run_cutadapt_umi(assembled_fastq, sample_dir, config, threads)
        if not umi_fastq:
            return None
        assembled_fastq = umi_fastq
        logger.info("Step 4: Running Bowtie2 and Samtools")
    else:
        logger.info("Step 3: Running Bowtie2 and Samtools")

    if not run_bowtie2_and_samtools(assembled_fastq, sample_dir, sample_name, config, threads):
        return None

    return sample_dir / f"{sample_name}.sorted.sam"


def cleanup_intermediate_files(sample_dir: Path, sample_name: Optional[str] = None) -> None:
    """
    Remove intermediate files generated during the pipeline for a sample,
    including the sorted SAM and merged FASTQ artifacts when they should not be kept.
    """
    files_to_remove = [
        sample_dir / 'fastp_R1.fastq.gz',
        sample_dir / 'fastp_R2.fastq.gz',
        sample_dir / 'pear.discarded.fastq.gz',
        sample_dir / 'pear.unassembled.forward.fastq.gz',
        sample_dir / 'pear.unassembled.reverse.fastq.gz',
        sample_dir / 'pear.trimmed.fastq',
        sample_dir / 'pear.trimmed.fastq.gz',
        sample_dir / 'pear.trimmed_with_flanks.fastq',
        sample_dir / 'pear.trimmed_with_flanks.fastq.gz',
    ]
    files_to_remove.extend(sample_dir.glob("*.assembled.fastq.gz"))
    if sample_name:
        files_to_remove.append(sample_dir / f"{sample_name}.sorted.sam")
    else:
        files_to_remove.extend(sample_dir.glob("*.sorted.sam"))

    for file_path in files_to_remove:
        try:
            if file_path.exists():
                file_path.unlink()
                logger.info(f"Removed intermediate file: {file_path}")
        except Exception as exc:
            logger.warning(f"Could not remove {file_path}: {exc}")


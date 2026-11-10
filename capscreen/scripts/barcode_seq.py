"""Barcode-seq assay path for CapScreen.

Shares FASTQ QC (FASTP/PEAR) with VP-seq, then aligns merged reads to a
Bowtie2 index built from a single barcode-amplicon FASTA (same pattern as
VP-seq: constant flanks with 15 Ns in the middle). Bowtie2 tolerates those
Ns. Counting extracts the DNA barcode between flanks (no peptide translation).

Build the index once, e.g.:
  bowtie2-build reference/barcodes.fasta reference/barcodes
Point barcode_seq.reference.genome / index_basename at that prefix.
"""

from __future__ import annotations

import copy
import logging
from pathlib import Path
from typing import Any, Dict, Optional

from capscreen.scripts import alignment as alignment_module
from capscreen.scripts import count as count_module

logger = logging.getLogger(__name__)

ASSAY_VP_SEQ = "vp-seq"
ASSAY_BARCODE_SEQ = "barcode-seq"
VALID_ASSAY_MODES = {ASSAY_VP_SEQ, ASSAY_BARCODE_SEQ}


def normalize_assay_mode(mode: Optional[str]) -> str:
    """Return a canonical assay mode string."""
    if mode is None or str(mode).strip() == "":
        return ASSAY_VP_SEQ
    value = str(mode).strip().lower().replace("_", "-")
    if value not in VALID_ASSAY_MODES:
        raise ValueError(
            f"Unknown assay mode '{mode}'. Expected one of: "
            f"{', '.join(sorted(VALID_ASSAY_MODES))}"
        )
    return value


def resolve_assay_config(config: Dict[str, Any], assay_mode: str) -> Dict[str, Any]:
    """
    Return a deep-copied config with flanking_sequences and reference set
    for the selected assay.

    VP-seq keeps top-level flanking_sequences/reference.
    Barcode-seq overlays barcode_seq.flanking_sequences and barcode_seq.reference.
    """
    mode = normalize_assay_mode(assay_mode)
    cfg = copy.deepcopy(config)
    cfg["assay_mode"] = mode

    if mode == ASSAY_VP_SEQ:
        return cfg

    barcode_cfg = cfg.get("barcode_seq") or {}
    flanks = barcode_cfg.get("flanking_sequences") or {}
    ref = barcode_cfg.get("reference") or {}

    if not flanks.get("flank_5p") or not flanks.get("flank_3p"):
        raise ValueError(
            "barcode-seq requires barcode_seq.flanking_sequences.flank_5p and flank_3p "
            "in config.json (constant sequence around the barcode NNNN insert)."
        )
    if not ref.get("genome"):
        raise ValueError(
            "barcode-seq requires barcode_seq.reference.genome "
            "(Bowtie2 -x index prefix built from barcodes.fasta, "
            "e.g. bowtie2-build barcodes.fasta barcodes)."
        )

    cfg["flanking_sequences"] = flanks
    merged_ref = dict(cfg.get("reference") or {})
    merged_ref.update(ref)
    cfg["reference"] = merged_ref
    return cfg


def default_library_path(config: Dict[str, Any], assay_mode: str) -> Optional[Path]:
    """Library CSV path from resolved assay config, if present."""
    cfg = resolve_assay_config(config, assay_mode)
    library = (cfg.get("reference") or {}).get("library")
    return Path(library) if library else None


def run_qc_and_alignment(
    fastq1: Path,
    fastq2: Path,
    output_dir: Path,
    sample_name: str,
    config: Dict[str, Any],
    threads: Optional[int],
    cut_umi: bool = False,
) -> Optional[Path]:
    """QC + align using barcode-seq flanks and barcode Bowtie2 reference."""
    logger.info("Assay mode: barcode-seq (align to barcode reference)")
    return alignment_module.run_qc_and_alignment(
        fastq1,
        fastq2,
        output_dir,
        sample_name,
        config,
        threads,
        cut_umi=cut_umi,
    )


def run_count(
    sam_file: Path,
    reference_file: Path,
    config: Dict[str, Any],
    output_file: Path,
    log_file: Optional[Path] = None,
    logger_instance: Optional[logging.Logger] = None,
    n_threads: Optional[int] = None,
) -> None:
    """Count DNA barcodes from SAM (no peptide translation)."""
    count_module.main(
        sam_file,
        reference_file,
        config,
        output_file,
        log_file=log_file,
        logger=logger_instance,
        n_threads=n_threads,
    )

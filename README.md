<p align="center">
  <img src="docs/capscreen-logo.jpg" alt="CapScreen, AAV capsid library screening" width="220">
</p>

<h1 align="center">CapScreen</h1>

<p align="center">
  AAV capsid library screening from FASTQ to pyDESeq2 and Er/Eh analysis.
</p>

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.18895265"><img src="https://zenodo.org/badge/1004327193.svg" alt="DOI"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="License: MIT"></a>
  <a href="https://github.com/hassansaei/CapScreen/releases"><img src="https://img.shields.io/badge/version-2.2.1-blue.svg" alt="Version 2.2.1"></a>
  <a href="https://hub.docker.com/r/saei/capscreen"><img src="https://img.shields.io/docker/v/saei/capscreen?label=docker&logo=docker" alt="Docker"></a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.9%E2%80%933.11-3776AB?logo=python&logoColor=white" alt="Python 3.9-3.11">
  <img src="https://img.shields.io/badge/Docker-ready-2496ED?logo=docker&logoColor=white" alt="Docker">
  <img src="https://img.shields.io/badge/stats-pyDESeq2-orange" alt="pyDESeq2">
  <img src="https://img.shields.io/badge/Apptainer-SIF-green" alt="Apptainer">
</p>

**CapScreen** is a pipeline designed for high-throughput screening of amplicon sequencing (NGS) results from AAV (Adeno-Associated Virus) library screens targeting specific receptors. The pipeline is suitable for both *binding* and *transduction* screening experiments, providing robust quality control, alignment, variant counting, and statistical analysis functionalities.


## Aim

The primary aim of CapScreen is to streamline the analysis of amplicon sequencing data generated from AAV library screens. It enables researchers to efficiently process, align, and quantify variants from large-scale binding or transduction screens, facilitating the identification of AAV variants with desired receptor targeting properties. The final output (matrix containing information of the reads, peptides and their counts) can be used for downstream statistical analysis and machine learning.


## Features

- **End-to-End Pipeline**: Automates all steps from raw FASTQ files to variant counts and statistical analysis.
- **Batch Processing**: Process multiple samples in parallel using a CSV file with sample metadata.
- **Quality Control**: Uses FASTP for adapter trimming, quality filtering, and reporting.
- **Read Merging**: Utilizes PEAR to merge paired-end reads, outputting compressed FASTQ files.
- **Alignment**: Aligns merged reads to a reference genome using Bowtie2 with customizable sensitivity and scoring.
- **Sorting and Indexing**: Employs Samtools for sorting and indexing alignment results.
- **Variant Extraction and Counting**: Extracts variable regions between user-defined flanking sequences, translates to peptides, and counts variant abundances.
- **Statistical Analysis**: Comprehensive statistical analysis including:
  - Differential expression with **pyDESeq2** (Python DESeq2; not the R/Bioconductor package)
  - RPM-based enrichment/enhancement (Er/Eh)
  - Principal Component Analysis (PCA)
  - Correlation analysis
  - Publication-ready plots (scatter plots, bar plots, heatmaps, PCA plots)
- **Flexible FASTQ Naming**: Supports multiple FASTQ naming conventions:
  - Illumina standard: `_1.fq.gz` / `_2.fq.gz` or `_1.fastq.gz` / `_2.fastq.gz`
  - Standard patterns: `_R1_` / `_R2_`, `_R1.` / `_R2.`, `_R1` / `_R2`
- **Comprehensive Logging**: Detailed logs for each step, with both console and file outputs.
- **Configurable**: All parameters are set via a single JSON configuration file.
- **Dockerized**: Ensures reproducibility and easy deployment across platforms.


## Workflow

The CapScreen pipeline consists of the following steps:

1. **Quality Control (FASTP)**
   - Trims adapters, filters low-quality reads, and generates QC reports.
2. **Read Merging (PEAR)**
   - Merges paired-end reads and outputs gzipped FASTQ files.
3. **Alignment (Bowtie2)**
   - Aligns merged reads to the reference genome.
4. **Sorting (Samtools)**
   - Sorts and indexes the alignment results.
5. **Variant Extraction and Counting**
   - Extracts the variable region between flanking sequences, translates to peptides, and counts each variant's abundance.
6. **Reporting**
   - Outputs a HTML report file, TSV file with variant counts and detailed logs.
7. **Statistical Analysis (Optional)**
   - Fits pairwise contrasts with pyDESeq2
   - Computes RPM-based enrichment/enhancement (Er/Eh)
   - Performs differential expression analysis
   - Generates PCA plots and correlation matrices
   - Creates publication-ready visualizations


## Usage

### 1. **Installation (Docker Recommended)**

CapScreen uses two Docker images:

| Image | Tag | Rebuild when |
| --- | --- | --- |
| Base (`ghcr.io/hassansaei/capscreen-base`) | `:latest` only (no version) | `conda/*.yml`, system tools, or `docker/DockerFile.base` change |
| App / final (`ghcr.io/hassansaei/capscreen`, and `saei/capscreen` if Hub secrets are set) | CapScreen version (e.g. `:2.2.1`) and `:latest` | Automatic on GitHub Release; or rebuild locally after code changes |

The base holds OS packages and the two Conda environments (`capscreen`, `capscreen-stat`).
The final image copies CapScreen source on top of that base. Publishing a GitHub Release
builds and pushes the app image via `.github/workflows/docker-app.yml`. Details:
[`docker/README.md`](docker/README.md).

**Rebuild the base** (infrequent; heavy):

```bash
docker build --platform linux/amd64 \
  -t ghcr.io/hassansaei/capscreen-base:latest \
  -f docker/DockerFile.base .
# then: docker push ghcr.io/hassansaei/capscreen-base:latest
```

**Rebuild the final app image** (usual path after code changes):

```bash
VERSION=$(python -c "import re; t=open('capscreen/version.py').read(); print(re.search(r\"__version__\\s*=\\s*['\\\"]([^'\\\"]+)['\\\"]\", t).group(1))")

docker build --platform linux/amd64 \
  -t saei/capscreen:${VERSION} \
  -t saei/capscreen:latest \
  -f docker/DockerFile .
```

Pull a published app image instead of building:

```bash
docker pull saei/capscreen:latest
```

### 2. **Assay modes (VP-seq vs barcode-seq)**

Choose the analysis path with a mutually exclusive flag (default: `--vp-seq`):

| Flag | QC | Alignment | Counting | Stats |
| --- | --- | --- | --- | --- |
| `--vp-seq` | FASTP + PEAR | Bowtie2 to VP/amplicon genome | Translate insert → peptides; match `ID,peptide` library | Same (`ID`) |
| `--barcode-seq` | FASTP + PEAR | Bowtie2 to **barcode FASTA index** | DNA barcode between flanks; match `ID,barcode` library | Same (`ID`) |

Barcode-seq uses `barcode_seq` in `config.json` for flanks and reference paths.
Provide **one** amplicon FASTA like VP-seq (`flank_5p` + **15 Ns** + `flank_3p`);
Bowtie2 tolerates the Ns. Then build indexes:

```bash
bowtie2-build reference/barcodes.fasta reference/barcodes
```

Point `barcode_seq.reference.genome` / `index_basename` at that prefix
(`barcodes.1.bt2`, …). Library CSV still lists real barcodes as
`ID,barcode`. Reference FASTA: `reference/barcodes.fasta`.
Library example: `reference/barcode_library.csv.example`.

```bash
capscreen pipeline --barcode-seq \
  --reference-file /opt/capscreen/reference/barcode_library.csv \
  ...
```

### 3. **Running the Pipeline**

#### Single Sample Mode

Mount your input/output directories and run the pipeline:

```bash
docker run --rm \
  -v /path/to/fastq:/opt/capscreen/input \
  -v /path/to/output:/opt/capscreen/output \
  capscreen \
  pipeline \
    --fastq1 /opt/capscreen/input/sample_R1.fastq.gz \
    --fastq2 /opt/capscreen/input/sample_R2.fastq.gz \
    --output-dir /opt/capscreen/output \
    --sample-name my_sample \
    --reference-file /opt/capscreen/input/reference_library.csv \
    --config /opt/capscreen/capscreen/config.json \
    --threads 8
```

#### Batch Mode (Multiple Samples)

Process multiple samples using a CSV file with sample metadata:

```bash
docker run --rm \
  -v /path/to/fastq:/opt/capscreen/input \
  -v /path/to/output:/opt/capscreen/output \
  capscreen \
  pipeline \
    --sample-info /opt/capscreen/input/Sample_info.csv \
    --output-dir /opt/capscreen/output \
    --reference-file /opt/capscreen/input/reference_library.csv \
    --merged-counts-output /opt/capscreen/output/merged.counts.csv \
    --config /opt/capscreen/capscreen/config.json \
    --threads 8 \
    --stat
```

**Sample_info.csv Format**

Required columns:
- `sample_id`: Unique identifier for each sample
- `run_mode`: Currently only `PE` (paired-end) is supported
- `path_FASTQ`: Path to R1 FASTQ file (R2 is automatically inferred)
- `group`: Sample group name (any label; e.g. `Target`, `MyKidney`, `Input`)
- `tech_rep`: Technical replicate number
- `bio_rep`: Biological replicate number
- `batch_seq`: Batch sequence number

Required for statistical analysis (`stat`: DESeq2 **and** enrichment/enhancement):
- `enrichment_role`: `input`, `target`, or `background` (blank = not used in DESeq2 or Er/Eh)
- `enrichment_arm`: Free-form arm tag that links one target group to its background group(s). A target used in several arms may list them separated by `;` (e.g. `OA_WT_vs_KO_POS;OA_WT_POS_vs_NEG`)

Group names may be anything. Contrasts are **not** inferred from CapScreen naming templates; they come only from these columns:
- Exactly one group must have `enrichment_role=input` (shared across arms; `enrichment_arm` may be blank)
- Each arm needs exactly one `target` group and one or more `background` groups with the same `enrichment_arm`
- All rows of the same `group` must share the same role and arm
- **DESeq2**: one contrast per target-versus-background pair (input is not a DESeq2 reference)
- **Enrichment/enhancement**: Er = target RPM / input RPM; Eh = target RPM / background RPM (all backgrounds for that arm pooled)
- Top-variant RPM bar plots show **all** groups in Sample_info
- PCA excludes the Sample_info input group(s)

Example for the sequencing pipeline only (roles optional until `stat` is run):
```csv
sample_id,run_mode,path_FASTQ,group,tech_rep,bio_rep,batch_seq
sample1,PE,/opt/capscreen/input/sample1_1.fq.gz,Target,1,1,1
sample2,PE,/opt/capscreen/input/sample2_1.fq.gz,Bead-Fc,1,1,1
sample3,PE,/opt/capscreen/input/sample3_1.fq.gz,Input,1,1,1
```

Example for `stat` (DESeq2 + enrichment/enhancement); roles are required:
```csv
sample_id,run_mode,path_FASTQ,group,tech_rep,bio_rep,batch_seq,enrichment_role,enrichment_arm
In1,PE,/opt/capscreen/input/In1_1.fq.gz,LibraryInput,1,1,1,input,
In2,PE,/opt/capscreen/input/In2_1.fq.gz,LibraryInput,1,2,1,input,
K1,PE,/opt/capscreen/input/K1_1.fq.gz,MyKidney,1,1,1,target,kidney_vs_offtarget
K2,PE,/opt/capscreen/input/K2_1.fq.gz,MyKidney,1,2,1,target,kidney_vs_offtarget
K3,PE,/opt/capscreen/input/K3_1.fq.gz,MyKidney,1,3,1,target,kidney_vs_offtarget
L1,PE,/opt/capscreen/input/L1_1.fq.gz,MyLiver,1,1,1,background,kidney_vs_offtarget
L2,PE,/opt/capscreen/input/L2_1.fq.gz,MyLiver,1,2,1,background,kidney_vs_offtarget
L3,PE,/opt/capscreen/input/L3_1.fq.gz,MyLiver,1,3,1,background,kidney_vs_offtarget
S1,PE,/opt/capscreen/input/S1_1.fq.gz,MySpleen,1,1,1,background,kidney_vs_offtarget
S2,PE,/opt/capscreen/input/S2_1.fq.gz,MySpleen,1,2,1,background,kidney_vs_offtarget
S3,PE,/opt/capscreen/input/S3_1.fq.gz,MySpleen,1,3,1,background,kidney_vs_offtarget
H1,PE,/opt/capscreen/input/H1_1.fq.gz,MyHeart,1,1,1,,
H2,PE,/opt/capscreen/input/H2_1.fq.gz,MyHeart,1,2,1,,
H3,PE,/opt/capscreen/input/H3_1.fq.gz,MyHeart,1,3,1,,
```

In that example, arm `kidney_vs_offtarget` drives:
- DESeq2 contrasts `MyKidney` vs `MyLiver` and `MyKidney` vs `MySpleen`
- Er = MyKidney / LibraryInput and Eh = MyKidney / (MyLiver + MySpleen)

`MyHeart` has a blank role, so it is skipped by DESeq2 and Er/Eh but still appears in the RPM bar plots.

**FASTQ Naming Conventions Supported:**
- Illumina standard: `sample_1.fq.gz` / `sample_2.fq.gz`
- Standard patterns: `sample_R1.fastq.gz` / `sample_R2.fastq.gz`
- Alternative: `sample_R1_001.fastq.gz` / `sample_R2_001.fastq.gz`

The tool automatically infers the R2 file from the R1 file path.

The `--stat` flag runs statistical analysis after creating the merged counts table.

The *reference_library.csv* is a user-specific file containing information about the variant IDs and peptide sequences of the input AAV library.

### 4. **Pipeline Subcommands**

- `pipeline`: Runs the full pipeline (QC, alignment, and counting). Supports both single sample and batch modes.
- `align`: Runs only QC and alignment steps.
- `count`: Runs only the variant counting step from a sorted SAM file.
- `report`: Generates only the HTML report for an existing sample directory.
- `stat`: Runs statistical analysis on merged counts table (requires merged.counts.csv and Sample_info.csv).

#### Examples

**Counting only:**
```bash
docker run --rm \
  -v /path/to/output:/data/output \
  -v /path/to/reference:/data/reference \
  capscreen \
  count \
    --sam-file /data/output/my_sample/my_sample.sorted.sam \
    --reference-file /data/reference/library.csv \
    --output /data/output/my_sample/sample.counts.csv \
    --output-dir /data/output \
    --sample-name my_sample \
    --config /opt/capscreen/capscreen/config.json
```

**Statistical analysis only:**
```bash
docker run --rm \
  -v /path/to/data:/opt/capscreen/input \
  -v /path/to/output:/opt/capscreen/output \
  capscreen \
  stat \
    --counts /opt/capscreen/input/merged.counts.csv \
    --sample-info /opt/capscreen/input/Sample_info.csv \
    --output-dir /opt/capscreen/output/statistics \
    --n-cpus 8 \
    --config /opt/capscreen/capscreen/config.json
```


## Configuration

All pipeline parameters are set in `capscreen/config.json`. Key sections include:

### Pipeline Configuration
- `fastp`: Quality control settings (compression, length requirements, duplicate calculation).
- `pear`: Read merging settings (overlap, assembly length, quality threshold).
- `bowtie2`: Alignment settings (mode, sensitivity, scoring parameters).
- `samtools`: Sorting and compression settings.
- `flanking_sequences`: 5' and 3' flanking sequences for variable region extraction.
- `reference`: Paths to reference genome and library.

Users can adapt `flanking_sequences` based on the AAV VP3 configuration and amplicon capture primers.

### Statistical Analysis Configuration

The `statistical_analysis` section configures downstream analysis:

- `group_roles`: Explicitly define group roles (optional, uses heuristics if not specified):
  - `input`: List of input group names
  - `target`: List of target group names
  - `background`: List of background group names

- `analysis`: Analysis parameters:
  - `n_cpus`: Number of CPUs for DESeq2
  - `lambda`: Pseudocount for RPM calculation
  - `low_expression_threshold`: Minimum expression threshold
  - `pca`: PCA parameters (components, random state, top features)

- `differential_expression`: DE analysis thresholds:
  - `padj_threshold`: Adjusted p-value threshold
  - `log2fc_threshold`: Log2 fold change threshold

- `output`: Result table formats:
  - `write_excel`: When true (default), dual-write DE and Er/Eh summary tables as TSV/CSV plus `.xlsx` companions (requires `openpyxl`)

- `plots`: Plotting configuration:
  - `scatter_highlight_top_n`: Number of top variants to highlight in scatter plots
  - `bar_plot_top_n`: Number of top variants in bar plots
  - `scatter.enable_highlighting`: Enable/disable highlighting in scatter plots (default: false)
  - `dpi`: Resolution settings for different plot types
  - `figure_sizes`: Figure dimensions for different plot types

- `colors`: Color scheme for plots (target, background, highlight colors) 

## Input/Output

### Input
- **Single Sample Mode**: 
  - Paired-end FASTQ files (R1 and R2)
  - Reference genome (Bowtie2 index)
  - Reference library CSV file (variant IDs and peptide sequences)
  
- **Batch Mode**:
  - Sample_info.csv file with sample metadata
  - Paired-end FASTQ files (R2 automatically inferred from R1)
  - Reference genome (Bowtie2 index)
  - Reference library CSV file

### Output
- **Per Sample**:
  - Processed and aligned files (SAM, sorted SAM)
  - Variant count table (TSV) with raw counts and RPM
  - HTML report with QC metrics and statistics
  - Detailed log files

- **Batch Mode**:
  - All per-sample outputs
  - Merged counts table (merged.counts.csv) with all samples
  - Sample_info.csv for statistical analysis

- **Statistical Analysis** (when `--stat` is used):
  - DESeq2 normalized counts
  - RPM enrichment analysis results
  - Differential expression results
  - PCA plots and correlation matrices
  - Presentation-ready plots (scatter, bar, heatmap, PCA)
  - Statistical analysis summary tables


## Best Practices

- **Use Docker**: Always use the provided Docker image for reproducibility and to avoid dependency issues.
- **Mount Directories**: Mount input/output directories as read/write volumes in Docker.
- **Configuration**: Review and customize `config.json` for your experiment, especially:
  - Flanking sequences (must match your amplicon design)
  - Bowtie2 sensitivity settings (adjust based on library diversity)
  - Statistical analysis parameters (group roles, thresholds)
- **Batch Processing**: Use batch mode (`--sample-info`) for processing multiple samples efficiently.
- **FASTQ Naming**: Ensure your FASTQ files follow supported naming conventions (see above).
- **Logging**: Check log files for detailed step-by-step information and troubleshooting.
- **Statistical Analysis**: 
  - Use `--stat` flag in batch mode to automatically run statistical analysis
  - Configure `group_roles` in config.json for predictable group classification
  - Review statistical outputs in the output directory
- **Resource Management**: 
  - Adjust `--threads` based on available CPU cores
  - Statistical analysis uses a separate Python 3.9–3.11 environment with pyDESeq2 (handled automatically in Docker)


## Citation

If you use CapScreen in your research, please cite it using the following DOI:

```bibtex
@software{saei_capscreen_2025,
  author       = {Saei, Hassan},
  title        = {{CapScreen: A pipeline for AAV library amplicon sequencing analysis}},
  year         = {2025},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18895265},
  url          = {https://doi.org/10.5281/zenodo.18895265}
}
```


## License

CapScreen is released under the [MIT License](LICENSE).

## Contact

For questions, issues, or contributions, please open an issue and provide detailed explanation of the issue possibly with the input file and running command. Thank you!

# CapScreen

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
  - DESeq2 normalization
  - RPM-based enrichment analysis
  - Differential expression analysis
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
   - Normalizes counts using DESeq2
   - Computes RPM-based enrichment
   - Performs differential expression analysis
   - Generates PCA plots and correlation matrices
   - Creates publication-ready visualizations


## Usage

### 1. **Installation (Docker Recommended)**

Build the Docker image from the project root:

```bash
docker build -f docker/DockerFile -t capscreen .
```

The Docker image includes two Conda environments:
- `capscreen`: Main pipeline environment (Python 3.7) for QC, alignment, and counting
- `capscreen-stat`: Statistical analysis environment (Python 3.9+) for DESeq2 and statistical analysis

### 2. **Running the Pipeline**

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
    --merged-counts-output /opt/capscreen/output/merged.counts.tsv \
    --config /opt/capscreen/capscreen/config.json \
    --threads 8 \
    --stat
```

**Sample_info.csv Format** (required columns):
- `sample_id`: Unique identifier for each sample
- `run_mode`: Currently only `PE` (paired-end) is supported
- `path_FASTQ`: Path to R1 FASTQ file (R2 is automatically inferred)
- `group`: Sample group name (e.g., "Target", "Bead-Fc", "Input")
- `tech_rep`: Technical replicate number
- `bio_rep`: Biological replicate number
- `batch_seq`: Batch sequence number

Example `Sample_info.csv`:
```csv
sample_id,run_mode,path_FASTQ,group,tech_rep,bio_rep,batch_seq
sample1,PE,/opt/capscreen/input/sample1_1.fq.gz,Target,1,1,1
sample2,PE,/opt/capscreen/input/sample2_1.fq.gz,Bead-Fc,1,1,1
sample3,PE,/opt/capscreen/input/sample3_1.fq.gz,Input,1,1,1
```

**FASTQ Naming Conventions Supported:**
- Illumina standard: `sample_1.fq.gz` / `sample_2.fq.gz`
- Standard patterns: `sample_R1.fastq.gz` / `sample_R2.fastq.gz`
- Alternative: `sample_R1_001.fastq.gz` / `sample_R2_001.fastq.gz`

The tool automatically infers the R2 file from the R1 file path.

The `--stat` flag runs statistical analysis after creating the merged counts table.

The *reference_library.csv* is a user-specific file containing information about the variant IDs and peptide sequences of the input AAV library.

### 3. **Pipeline Subcommands**

- `pipeline`: Runs the full pipeline (QC, alignment, and counting). Supports both single sample and batch modes.
- `align`: Runs only QC and alignment steps.
- `count`: Runs only the variant counting step from a sorted SAM file.
- `report`: Generates only the HTML report for an existing sample directory.
- `stat`: Runs statistical analysis on merged counts table (requires merged.counts.tsv and Sample_info.csv).

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
    --output /data/output/my_sample/counts.tsv \
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
    --counts /opt/capscreen/input/merged.counts.tsv \
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
  - Merged counts table (merged.counts.tsv) with all samples
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
  - Statistical analysis uses a separate Python 3.9+ environment (handled automatically in Docker)


## Citation

Coming soon....


## License

MIT License

## Contact

For questions, issues, or contributions, please open an issue and provide detailed explanation of the issue possibly with the input file and running command. Thank you!

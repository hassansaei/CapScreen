# CapScreen

**CapScreen** is a pipeline designed for high-throughput screening of amplicon sequencing (NGS) results from AAV (Adeno-Associated Virus) library screens targeting specific receptors. The pipeline is suitable for both *binding* and *transduction* screening experiments, providing robust quality control, alignment, and variant counting functionalities. More options will be added overtime.


## Aim

The primary aim of CapScreen is to streamline the analysis of amplicon sequencing data generated from AAV library screens. It enables researchers to efficiently process, align, and quantify variants from large-scale binding or transduction screens, facilitating the identification of AAV variants with desired receptor targeting properties. The final output (matrix containing inforation of the reads. peptides and their count) will be used for downstream ML analysis.


## Features

- **End-to-End Pipeline**: Automates all steps from raw FASTQ files to variant counts.
- **Quality Control**: Uses FASTP for adapter trimming, quality filtering, and reporting.
- **Read Merging**: Utilizes PEAR to merge paired-end reads, outputting compressed FASTQ files.
- **Alignment**: Aligns merged reads to a reference genome using Bowtie2 with customizable sensitivity and scoring.
- **Sorting and Indexing**: Employs Samtools for sorting and indexing alignment results.
- **Variant Extraction and Counting**: Extracts variable regions between user-defined flanking sequences, translates to peptides, and counts variant abundances.
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
   - Outputs a TSV file with variant counts and detailed logs.


## Usage

### 1. **Installation (Docker Recommended)**

Build the Docker image from the project root:

```bash
docker build -f docker/DockerFile -t capscreen .
```

### 2. **Running the Pipeline**

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

The *reference_library.csv* is a user-specific file containing information about the variant IDs and peptide sequences of the input AAV library. 


### 3. **Pipeline Subcommands**

- `pipeline`: Runs the full pipeline (QC, alignment, and counting).
- `align`: Runs only QC and alignment steps.
- `count`: Runs only the variant counting step from a sorted SAM file.

Example for counting only:

```bash
docker run --rm \
  -v /path/to/output:/data/output \
  -v /path/to/reference:/data/reference \
  capscreen \
  count \
    /data/output/my_sample/aligned.sorted.sam \
    /data/reference/library.csv \
    --output /data/output/my_sample/counts.tsv \
    --config /opt/capscreen/capscreen/config.json
```


## Configuration

All pipeline parameters are set in `capscreen/config.json`. Key sections include:

- `fastp`: Quality control settings.
- `pear`: Read merging settings.
- `bowtie2`: Alignment settings.
- `samtools`: Sorting and compression settings.
- `flanking_sequences`: 5' and 3' flanking sequences for variable region extraction.
- `reference`: Paths to reference genome and library.

User can adapt falnking_sequence based on the AAV VP3 configuration and amplicon capture primers. 

## Input/Output

- **Input**: Paired-end FASTQ files, reference genome (Bowtie2 index), and reference library (CSV).
- **Output**: Processed and aligned files, variant count tables (TSV), and log files.


## Best Practices

- Always use the provided Docker image for reproducibility.
- Mount input/output directories as read/write volumes.
- Review and customize `config.json` for your experiment.
- Check log files for detailed step-by-step information and troubleshooting.


## Citation

Coming soon....


## License

MIT License

## Contact

For questions, issues, or contributions, please open an issue and provide detailed explanation of the issue possibly with the input file and running command. Thank you!

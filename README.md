Of course. Here is the updated `README.md`.

This version has been thoroughly revised to reflect the latest pipeline capabilities, including:

  * Handling of both **paired-end and single-end** data.
  * The new, more detailed **samplesheet format**.
  * The updated **output structure**, which saves key results directly alongside the source data.
  * A more accurate description of the current **workflow logic**.

All typos like `mani.nf` have been corrected to `main.nf`.

-----

# RNA-Seq Analysis Pipeline

\<p align="center"\>
\<a href="[suspicious link removed]" target="\_blank" rel="noopener noreferrer"\>\<img src="[suspicious link removed]" alt="Nextflow"\>\</a\>
\<a href="[suspicious link removed]" target="\_blank" rel="noopener noreferrer"\>\<img src="[suspicious link removed]" alt="License"\>\</a\>
\</p\>

A robust and flexible Nextflow pipeline for the analysis of RNA-sequencing data. It supports both paired-end and single-end reads, performing quality control, adapter trimming, alignment, and gene-level quantification.

-----

## Table of Contents

  - [Introduction](https://www.google.com/search?q=%23introduction)
  - [Pipeline Workflow](https://www.google.com/search?q=%23pipeline-workflow)
  - [System Requirements](https://www.google.com/search?q=%23system-requirements)
  - [Installation & Setup](https://www.google.com/search?q=%23installation--setup)
  - [Usage](https://www.google.com/search?q=%23usage)
  - [Output Directory Structure](https://www.google.com/search?q=%23output-directory-structure)
  - [Troubleshooting](https://www.google.com/search?q=%23troubleshooting)
  - [Contributing](https://www.google.com/search?q=%23contributing)
  - [License](https://www.google.com/search?q=%23license)
  - [Contact](https://www.google.com/search?q=%23contact)

-----

## Introduction

This pipeline processes raw FASTQ files through a series of standard bioinformatics tools to produce a gene count matrix and a comprehensive quality control report. Built with Nextflow, it offers excellent scalability and reproducibility, enabling it to run on local machines, HPC clusters, or cloud environments with minimal configuration changes.

### Key Features

  - **Handles Both SE & PE Data**: Processes paired-end and single-end data seamlessly based on a simple samplesheet annotation.
  - **Flexible Execution Modes**:
    1.  **Full Pipeline Mode**: Executes the complete workflow from trimming to counting.
    2.  **Trim-Only Mode**: Runs only trimming and MultiQC to assess trimming results before a full run.
    3.  **QC-Only Mode**: Runs an initial `FastQC` and `MultiQC` on raw reads to assess data quality.
  - **Targeted Outputs**: Saves trimmed reads and alignments in a `trimmed/` subdirectory alongside the original data for easy access.
  - **Portability**: Supports various execution environments like SLURM, Conda, Docker, and Singularity.
  - **Reproducibility**: Ensures that the analysis is reproducible by managing dependencies and workflow versions.
  - **Customizable**: Per-sample trimming parameters can be specified directly in the samplesheet.

-----

## Pipeline Workflow

The pipeline consists of the following major steps:

1.  **Input Reading**: The `samplesheet.csv` is parsed to identify sample metadata and read paths.
2.  **Data Branching**: Samples are automatically branched into **single-end** or **paired-end** workflows based on the `trim_type` column.
3.  **SE Read Combination (CAT\_SE\_READS)**: For single-end samples with multiple FASTQ files (e.g., from different lanes like L6/L7), the files are concatenated into one.
4.  **Adapter/Quality Trimming (Trim Galore)**: Adapters and low-quality bases are removed from each sample.
5.  **Quality Control (FastQC)**: Per-sample quality control reports are generated from the trimmed reads.
6.  **Genome Indexing (HISAT2)**: A genome index is built from a reference FASTA file (this step is skipped on subsequent runs).
7.  **Alignment (HISAT2)**: Trimmed reads (both SE and PE) are aligned to the reference genome.
8.  **Gene Quantification (featureCounts)**: Aligned reads are assigned to genes based on a GTF annotation file. PE and SE samples are processed separately to ensure correct counting parameters.
9.  **Merge Counts (MERGE\_COUNTS)**: The count results from the PE and SE streams are merged into a final gene count matrix.
10. **Aggregate Reporting (MultiQC)**: Results and logs from all tools are combined into a single, interactive HTML report.

-----

## System Requirements

  - **Nextflow**: Version `21.10.3` or higher.
  - **Execution Environment**:
      - An HPC system with a job scheduler like SLURM.
      - Alternatively, Conda, Docker, or Singularity can be used for dependency management.
  - **Dependencies**:
      - Trim Galore
      - FastQC
      - HISAT2
      - SAMtools
      - Subread/featureCounts
      - Python (for `MULTIQC` and helper scripts)

-----

## Installation & Setup

1.  **Clone the Repository**:

    ```bash
    git clone <repository_url>
    cd <repository_name>
    ```

2.  **Configure the Pipeline**:
    Edit the `nextflow.config` file to specify paths for your reference genome and annotation files.

    ```groovy
    params {
        // ... other params
        genome_fasta = '/path/to/your/genome/hg38.fa'
        genome_gtf = '/path/to/your/genome/hg38.gtf'
        // ...
    }
    ```

-----

## Usage

### 1\. Prepare Input Data

Create a `samplesheet.csv` file with the structure detailed below. This sheet is critical for defining samples, their read types, and any custom parameters.

**`samplesheet.csv` Columns:**

| Column | Description | Example |
| :--- | :--- | :--- |
| `sample` | New, clean sample identifier. Used for naming output files. | `Hcy_293T_Rep1` |
| `orig_name_1` | Original filename of the first read file for traceability. | `Hcy_1_293T_4_1.fq.gz` |
| `orig_name_2` | Original filename of the second read file. | `Hcy_1_293T_4_2.fq.gz` |
| `fastq_1` | Full path to the first FASTQ file (Read 1 for PE, or first SE file). | `/path/to/Hcy_1_293T_4_1.fq.gz` |
| `fastq_2` | Full path to the second FASTQ file (Read 2 for PE, or second SE file). | `/path/to/Hcy_1_293T_4_2.fq.gz` |
| `condition` | Experimental condition metadata. | `Hcy` |
| `cell_line`| Cell line metadata. | `293T` |
| `parental_line`| Parental cell line metadata. | `293T` |
| `trim_type`| **Crucial**: Set to `paired` or `single` to direct the workflow. | `paired` |
| `trim_args`| Optional arguments to override default Trim Galore settings. | `--quality 15` |

**Example `samplesheet.csv`**:

```csv
sample,orig_name_1,orig_name_2,fastq_1,fastq_2,condition,cell_line,parental_line,trim_type,trim_args
Hcy_293T_Rep1,Hcy_1_293T_4_1.fq.gz,Hcy_1_293T_4_2.fq.gz,/path/to/data/Hcy_1_293T_4_1.fq.gz,/path/to/data/Hcy_1_293T_4_2.fq.gz,Hcy,293T,293T,paired,
Hcy_468_Rep1,4R044-L6-P19.gz,4R044-L7-P19.gz,/path/to/data/Hcy_1_468_12hr_19/4R044-L6-P19.gz,/path/to/data/Hcy_1_468_12hr_19/4R044-L7-P19.gz,Hcy,468,468,single,--quality 15 --clip_R1 5
```

### 2\. Run the Pipeline

Navigate to the pipeline directory and execute one of the following commands. The main script is **`main.nf`**. We recommend using the `hpc` profile for cluster execution.

#### Mode 1: Full Pipeline Run

Runs the complete workflow from trimming to feature counting.

```bash
nextflow run main.nf -profile hpc --outdir ./results
```

#### Mode 2: Trim-Only Run

Runs only `Trim Galore` and `MultiQC` to check trimming results.

```bash
nextflow run main.nf -profile hpc --trim_only true --outdir ./results
```

#### Mode 3: QC-Only Run

Runs `FastQC` and `MultiQC` on the **raw reads** for an initial quality assessment.

```bash
nextflow run main.nf -profile hpc --qc_only true --outdir ./results
```

### Resuming the Pipeline

If the pipeline is interrupted, you can resume it from the last successfully completed step using the `-resume` flag.

```bash
nextflow run main.nf -profile hpc --outdir ./results -resume
```

### Command-line Parameters

| Parameter | Description |
| :--- | :--- |
| `--outdir` | The directory where central results (counts, reports) will be saved (Default: `results`). |
| `--samplesheet_file` | Path to the input samplesheet (Default: `samplesheet.csv`). |
| `--trim_only` | A boolean flag to activate Trim-only mode (Default: `false`). |
| `--qc_only` | A boolean flag to activate QC-only mode on raw reads (Default: `false`). |

-----

## Output Directory Structure

The pipeline produces outputs in two locations: a central results folder and a `trimmed` subfolder alongside your original data.

### 1\. Per-Sample Outputs (In Source Data Directory)

For each sample, a new directory named `trimmed/` is created within its original folder. This makes it easy to find processed files related to the source data.

```
/path/to/your/data/
└── Hcy_1_293T_4/
    ├── Hcy_1_293T_4_1.fq.gz          (Original Input)
    ├── Hcy_1_293T_4_2.fq.gz          (Original Input)
    └── trimmed/                      (NEW OUTPUT FOLDER)
        ├── Hcy_293T_Rep1.bam
        ├── Hcy_293T_Rep1.bam.bai
        ├── Hcy_293T_Rep1.summary.txt
        ├── Hcy_293T_Rep1_trimming_report.txt
        ├── Hcy_293T_Rep1_val_1.fq.gz
        └── Hcy_293T_Rep1_val_2.fq.gz
```

### 2\. Central Output Directory (Defined by `--outdir`)

This directory contains pipeline-wide results and reports.

```
results/
├── counts/
│   ├── counts.tsv                  # Final merged gene count matrix
│   └── counts.tsv.summary          # Summary of featureCounts run
├── fastqc_raw/                       # FastQC reports (QC-only mode)
│   └── ...
├── genome_index/
│   └── ...
└── multiqc_report.html               # Final aggregated MultiQC report
```

-----

## Troubleshooting

  - **Issue**: Pipeline fails with an error related to environment modules.

      - **Solution**: Ensure that the required software is available as environment modules on your HPC system and that the module names in `nextflow.config` are correct.

  - **Issue**: `featureCounts` fails due to mixed SE/PE bams.

      - **Solution**: This pipeline is designed to prevent this by separating SE and PE samples into different `featureCounts` processes before merging the text-based results. If this error occurs, check the `trim_type` column in your samplesheet for accuracy.

-----

## Contributing

Contributions are welcome\! Please feel free to submit a pull request or open an issue to report bugs or suggest improvements.

-----

## License

This project is licensed under the MIT License.

-----

## Contact

For questions or support, please contact David at [tealeave@gmail.com].
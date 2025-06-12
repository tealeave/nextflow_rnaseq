Of course. Based on the robust and flexible Nextflow script we've built, here is a comprehensive `README.md` file. You can save this directly into your project's directory.

---

# Flexible RNA-Seq Analysis Pipeline

This repository contains a modular and resource-aware RNA-Seq analysis pipeline built with Nextflow. It is designed for use in HPC environments (specifically with the SLURM scheduler) and provides granular control over execution steps to manage storage and computational resources effectively.

The pipeline performs standard RNA-Seq analysis, including quality control, adapter trimming, alignment, and gene-level quantification.

## Key Features

-   **Step-wise Execution**: Run the entire pipeline or stop after a specific stage (QC, trimming, or alignment).
-   **Storage Management**: Optimized for environments with limited home directory space by allowing the large temporary `work` directory to be placed on a separate, larger storage volume.
-   **Efficient File Handling**: Uses `move` for publishing results, cleaning up the work directory as it runs to conserve space.
-   **Flexible Input**: Handles both single-end and paired-end data, as well as samples with multiple FASTQ files (e.g., from different lanes).
-   **Customizable**: Per-sample trimming parameters can be specified directly in the samplesheet.

## Pipeline Steps

The pipeline uses the following tools:

1.  **FastQC** (`v0.11.9`): Initial quality control on raw and trimmed reads.
2.  **Trim Galore!** (`v0.6.7`): Adapter and quality trimming.
3.  **HISAT2** (`v2.2.1`): Alignment of reads to a reference genome.
4.  **Samtools** (`v1.15.1`): Sorting and indexing of BAM alignment files.
5.  **featureCounts** (Subread `v2.0.3`): Gene-level quantification.
6.  **MultiQC**: Aggregates analysis results and logs into a final, unified report.

## Prerequisites

Before running the pipeline, please ensure you have the following installed and available in your environment:

-   **Nextflow** (`~21.10.x` or later)
-   **HPC Environment Modules**: The pipeline is configured to use environment modules for its software dependencies. Ensure that modules for the tools listed above are available on your system.

## Setup

1.  **Clone or download the pipeline files:**
    Place `main.nf` and `nextflow.config` in your project directory.

2.  **Prepare Reference Genomes:**
    You will need a reference genome in FASTA format and a corresponding gene annotation file in GTF format. You can specify their paths during pipeline execution.

3.  **Create a Samplesheet:**
    The pipeline requires a samplesheet in CSV format (`.csv`) that details the input files. The default filename is `samplesheet.csv`.

    **Columns:**
    | Column | Description | Required |
    | :--- | :--- | :--- |
    | `sample` | A unique identifier for the sample. No spaces or special characters. | **Yes** |
    | `fastq_1` | Full path to the forward read file (`_R1.fastq.gz`). For single-end data with multiple files per sample, this should be the first file. | **Yes** |
    | `fastq_2` | Full path to the reverse read file (`_R2.fastq.gz`). **Leave this column empty for single-end data.** | No |
    | `trim_type` | The type of sequencing data. Must be either `single` or `paired`. | **Yes** |
    | `trim_args`| Optional custom arguments for Trim Galore!. If left blank, the default from `nextflow.config` is used. | No |

    **Example `samplesheet.csv`:**
    ```csv
    sample,fastq_1,fastq_2,trim_type,trim_args
    Sample_A_PE,/path/to/data/sampA_R1.fq.gz,/path/to/data/sampA_R2.fq.gz,paired,
    Sample_B_PE,/path/to/data/sampB_R1.fq.gz,/path/to/data/sampB_R2.fq.gz,paired,--quality 15 --length 25
    Sample_C_SE,/path/to/data/sampC.fq.gz,,single,
    ```

## Usage

### Basic Command

The basic command to execute the pipeline is:

```bash
nextflow run main.nf -profile hpc [options]
```

### Core Parameters

| Parameter | Description | Example |
| :--- | :--- | :--- |
| `--samplesheet_file` | Path to the input samplesheet CSV file. | `--samplesheet_file samples.csv` |
| `--outdir` | Path to the directory where results will be saved. | `--outdir ./results` |
| `--genome_fasta` | Path to the reference genome FASTA file. | `--genome_fasta /refs/hg38.fa` |
| `--genome_gtf` | Path to the genome annotation GTF file. | `--genome_gtf /refs/hg38.gtf` |

### Storage Management

To prevent filling up your local or home storage, you can specify a location for the temporary `work` directory on a larger storage volume using the `-w` flag.

```bash
# Example: Use a work directory on a large shared drive
nextflow run main.nf -profile hpc -w /path/to/large/storage/work_dir
```

### Step-wise Execution

Use the `--step` parameter to control how far the pipeline runs. This is useful for debugging, resource management, and milestone checks.

-   `--step qc_raw`: Generates FastQC reports on raw reads and stops.
-   `--step trim`: Runs `qc_raw` steps, then trims reads, generates FastQC reports on trimmed reads, and stops.
-   `--step align`: Runs `trim` steps, then aligns reads to the genome, and stops.
-   `--step full` (Default): Runs the complete pipeline, including feature counting.

### Resuming a Pipeline

If the pipeline is interrupted, you can resume it from the last successful step using the `-resume` flag. Nextflow will use cached results for completed tasks.

```bash
# If the pipeline stopped during the alignment step, resume it
nextflow run main.nf -profile hpc -w /path/to/work_dir -resume
```

### Example Commands

**1. Full Run on HPC**
Run the complete pipeline using the `hpc` profile and specify a custom work directory and output folder.

```bash
nextflow run main.nf -profile hpc \
    -w /path/to/large_storage/nextflow_work \
    --samplesheet_file samplesheet.csv \
    --outdir ./final_results \
    --genome_fasta /refs/hg38.fa \
    --genome_gtf /refs/hg38.gtf
```

**2. Run Trimming Step Only**
Execute the pipeline only up to the trimming and post-trimming QC stage.

```bash
nextflow run main.nf -profile hpc -w /path/to/work_dir --step trim
```

**3. Resume and Run Alignment**
If you have already completed the `trim` step, you can resume the pipeline to run the `align` step.

```bash
nextflow run main.nf -profile hpc -w /path/to/work_dir --step align -resume
```

## Output Directory Structure

The pipeline will create an output directory (default: `results/`) with the following structure:

```
results/
├── alignment/
│   ├── Sample_A_PE/
│   │   ├── Sample_A_PE.bam
│   │   ├── Sample_A_PE.bam.bai
│   │   └── Sample_A_PE.summary.txt
│   └── ...
├── counts/
│   ├── counts.tsv
│   └── counts.tsv.summary
├── fastqc_raw/
│   ├── ...
├── fastqc_trimmed/
│   ├── ...
├── genome_index/
│   ├── ...
├── trimmed_reads/
│   ├── Sample_A_PE/
│   │   ├── Sample_A_PE_val_1.fq.gz
│   │   ├── Sample_A_PE_val_2.fq.gz
│   │   └── Sample_A_PE_trimming_report.txt
│   └── ...
└── multiqc_report.html
```
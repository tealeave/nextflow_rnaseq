#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- Main Workflow ---
workflow {

    // --- Input Channel ---
    // Reads the samplesheet.csv, creates a map for each row,
    // and prepares the [meta, [fastq_1, fastq_2]] tuple for downstream processes.
    Channel
        .fromPath(params.samplesheet_file)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def fastq_1 = file(row.fastq_1)
            def fastq_2 = file(row.fastq_2)
            if (!fastq_1.exists() || !fastq_2.exists()) {
                error "FASTQ file not found for sample ${meta.id}:\n  ${row.fastq_1}\n  ${row.fastq_2}"
            }
            return [meta, [fastq_1, fastq_2]]
        }
        .ifEmpty { error "Samplesheet is empty or not found: ${params.samplesheet_file}" }
        .set { ch_reads }

    // --- Step 1: Quality Control (FastQC) ---
    ch_fastqc_results = FASTQC(ch_reads)

    // --- Step 2: Build HISAT2 Genome Index ---
    // This process runs only once.
    ch_hisat2_index = HISAT2_INDEX(file(params.genome_fasta))

    // --- Step 3: Alignment (HISAT2) ---
    // Aligns reads for each sample to the reference genome.
    ch_bams = HISAT2_ALIGN(ch_reads, ch_hisat2_index)

    // --- Step 4: Generate Count Matrix (featureCounts) ---
    // Collects all BAM files and creates a single gene count matrix.
    ch_featurecounts = FEATURE_COUNTS(ch_bams.map{ it[1] }.collect(), file(params.genome_gtf))

    // --- Step 5: Aggregate All Results (MultiQC) ---
    // Collects logs and reports from all previous steps into a single HTML report.
    ch_multiqc = MULTIQC(
        ch_fastqc_results.collect(),
        ch_bams.map{ it[2] }.collect(), // Alignment summaries
        ch_featurecounts.summary.collect()
    )
}

// ========================================================================================
//                              PROCESS DEFINITIONS
// ========================================================================================

process FASTQC {
    tag "$meta.id"
    publishDir "${params.outdir}/fastqc", mode: 'copy', pattern: '*.{html,zip}'

    input:
    tuple val(meta), path(reads)

    output:
    path("*.{html,zip}"), emit: report

    script:
    """
    fastqc -o . -q ${reads.join(' ')}
    """
}

process HISAT2_INDEX {
    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    path genome_fasta

    output:
    path("hisat2_index")

    script:
    """
    mkdir hisat2_index
    hisat2-build -p ${task.cpus} ${genome_fasta} hisat2_index/genome
    """
}

process HISAT2_ALIGN {
    tag "$meta.id"
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*.{bam,bai,summary.txt}"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path("*.summary.txt"), emit: summary
    path("*.bai"), emit: bai // for downstream convenience

    script:
    def prefix = "${meta.id}"
    """
    hisat2 -p ${task.cpus} \\
        --dta \\
        --rna-strandness RF \\
        --summary-file ${prefix}.summary.txt \\
        -x ${index}/genome \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} | \\
    samtools view -bS - > ${prefix}.unsorted.bam

    samtools sort -@ ${task.cpus} -o ${prefix}.bam ${prefix}.unsorted.bam
    samtools index ${prefix}.bam
    """
}

process FEATURE_COUNTS {
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path bams
    path gtf

    output:
    path "counts.tsv", emit: main
    path "counts.tsv.summary", emit: summary

    script:
    """
    featureCounts -p -s 2 -T ${task.cpus} \\
        -a ${gtf} \\
        -o counts.tsv \\
        ${bams.join(' ')}
    """
}

process MULTIQC {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}
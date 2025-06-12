#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ========================================================================================
//                                  P I P E L I N E  S U M M A R Y
// ========================================================================================
// This pipeline performs RNA-seq analysis with three modes:
// 1. Full Mode (default): Concatenates SE reads, trims all reads, runs FastQC, HISAT2
//    alignment, featureCounts, and MultiQC.
// 2. QC-Only Mode: Runs only FastQC on raw reads and MultiQC. Activate with `--qc_only true`.
// 3. Trim-Only Mode: Combines and trims reads, then runs MultiQC. Activate with `--trim_only true`.
// ========================================================================================


// --- Main Workflow ---
workflow {
    // --- Input Channel ---
    Channel
        .fromPath(params.samplesheet_file)
        .splitCsv(header: true)
        .map { row ->
            def meta = row.clone()
            meta.id = row.sample
            meta.orig_dir = file(row.fastq_1).getParent() // Store original directory for publishing results

            if (!meta.trim_args?.trim()) {
                meta.trim_args = params.trim_args_default
            }

            def file_list = [ file(row.fastq_1) ]
            if (row.fastq_2?.trim()) { // Only add fastq_2 if it's specified
                file_list.add(file(row.fastq_2))
            }
            return [meta, file_list]
        }
        .ifEmpty { error "Samplesheet is empty or not found: ${params.samplesheet_file}" }
        .set { ch_raw_reads }

    // --- Branch raw reads based on trim_type for processing ---
    ch_raw_reads
        .branch { meta, files ->
            single: meta.trim_type == 'single'
            paired: meta.trim_type == 'paired'
            other:  true
        }
        .set { ch_branched_reads }

    ch_branched_reads.other
        .count()
        .subscribe { count ->
            if (count > 0) {
                error "${count} samples found with 'trim_type' other than 'single' or 'paired'. Please check samplesheet."
            }
        }

    // --- Process Single-End Reads ---
    CAT_SE_READS(ch_branched_reads.single)
    TRIM_GALORE_SE(CAT_SE_READS.out.reads)

    // --- Process Paired-End Reads ---
    TRIM_GALORE_PE(ch_branched_reads.paired)

    // --- Combine outputs for downstream steps ---
    ch_trimmed_reads = TRIM_GALORE_SE.out.trimmed_reads.mix(TRIM_GALORE_PE.out.trimmed_reads)
    ch_trim_reports  = TRIM_GALORE_SE.out.report.mix(TRIM_GALORE_PE.out.report)

    // --- Conditional Workflow Execution ---
    if (params.trim_only) {
        log.info "Running in Trim-only mode."
        MULTIQC(ch_trim_reports.collect().ifEmpty { [] })

    } else if (params.qc_only) {
        log.info "Running in QC-only mode (on raw reads)."
        FASTQC(ch_raw_reads)
        MULTIQC(FASTQC.out.report.collect().ifEmpty { [] })

    } else {
        // --- Full Pipeline Mode ---
        log.info "Running the full pipeline."

        FASTQC(ch_trimmed_reads)
        HISAT2_INDEX(file(params.genome_fasta))
        HISAT2_ALIGN(ch_trimmed_reads, HISAT2_INDEX.out)

        // Branch BAM files for correct feature counting
        HISAT2_ALIGN.out.bam
            .branch { meta, bam ->
                single: meta.trim_type == 'single'
                paired: meta.trim_type == 'paired'
            }
            .set { ch_bams }

        FEATURE_COUNTS_PE(ch_bams.paired.map{ it[1] }.collect().ifEmpty{[]}, file(params.genome_gtf))
        FEATURE_COUNTS_SE(ch_bams.single.map{ it[1] }.collect().ifEmpty{[]}, file(params.genome_gtf))

        MERGE_COUNTS(
            FEATURE_COUNTS_PE.out.main.ifEmpty(file('empty.tsv')),
            FEATURE_COUNTS_PE.out.summary.ifEmpty(file('empty.summary')),
            FEATURE_COUNTS_SE.out.main.ifEmpty(file('empty.tsv')),
            FEATURE_COUNTS_SE.out.summary.ifEmpty(file('empty.summary'))
        )

        def multiqc_inputs = Channel.empty()
            .mix(ch_trim_reports)
            .mix(FASTQC.out.report)
            .mix(HISAT2_ALIGN.out.summary)
            .mix(MERGE_COUNTS.out.summary)
            .collect()
        MULTIQC(multiqc_inputs)
    }
}

// ========================================================================================
//                                P R O C E S S  D E F I N I T I O N S
// ========================================================================================

process FASTQC {
    tag "$meta.id"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy', pattern: '*.{html,zip}'

    input:
    tuple val(meta), path(reads)

    output:
    path("*.{html,zip}"), emit: report

    script:
    """
    fastqc -o . -q ${reads.join(' ')}
    """
}

process CAT_SE_READS {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads) // Expecting [L6.fq.gz, L7.fq.gz]

    output:
    tuple val(meta), path("${meta.id}.fq.gz"), emit: reads

    script:
    """
    cat ${reads.join(' ')} > ${meta.id}.fq.gz
    """
}

process TRIM_GALORE_PE {
    tag "$meta.id"
    publishDir path: "${meta.orig_dir}/trimmed", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_val_*.fq.gz"), emit: trimmed_reads
    path("${meta.id}_trimming_report.txt"),        emit: report

    script:
    """
    trim_galore --paired \\
        --basename ${meta.id} \\
        --cores ${task.cpus} \\
        ${meta.trim_args} \\
        --output_dir . \\
        ${reads[0]} ${reads[1]}
    """
}

process TRIM_GALORE_SE {
    tag "$meta.id"
    publishDir path: "${meta.orig_dir}/trimmed", mode: 'copy'

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fq.gz"), emit: trimmed_reads
    path("${meta.id}_trimming_report.txt"),          emit: report

    script:
    def se_trim_args = meta.trim_args.replaceAll('--clip_R2 \\d+', '').trim()
    """
    trim_galore \\
        --basename ${meta.id} \\
        --cores ${task.cpus} \\
        ${se_trim_args} \\
        --output_dir . \\
        ${read}
    """
}

process HISAT2_INDEX {
    publishDir "${params.outdir}/genome_index", mode: 'copy'

    input:
    path genome_fasta

    output:
    path("hisat2_index"), emit: index

    script:
    def prefix = genome_fasta.baseName
    """
    mkdir hisat2_index
    hisat2-build -p ${task.cpus} ${genome_fasta} hisat2_index/${prefix}
    """
}

process HISAT2_ALIGN {
    tag "$meta.id"
    publishDir path: "${meta.orig_dir}/trimmed", mode: 'copy', pattern: "${meta.id}*"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    path("${meta.id}.summary.txt"),          emit: summary
    path("${meta.id}.bam.bai"),              emit: bai

    script:
    def index_prefix = file(index).list('*.ht2')[0].baseName
    def read_inputs = (meta.trim_type == 'paired') ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads[0]}"
    """
    hisat2 -p ${task.cpus} \\
        --dta \\
        --rna-strandness RF \\
        --summary-file ${meta.id}.summary.txt \\
        -x ${index}/${index_prefix} \\
        ${read_inputs} | \\
    samtools view -bS - > ${meta.id}.unsorted.bam

    samtools sort -@ ${task.cpus} -o ${meta.id}.bam ${meta.id}.unsorted.bam
    samtools index ${meta.id}.bam
    """
}

process FEATURE_COUNTS_PE {
    publishDir "${params.outdir}/counts", mode: 'copy', pattern: "counts_pe.tsv*"

    input:
    path bams
    path gtf

    output:
    path "counts_pe.tsv",         emit: main
    path "counts_pe.tsv.summary", emit: summary

    when:
    bams

    script:
    """
    featureCounts -p -s 2 -T ${task.cpus} \\
        -a ${gtf} \\
        -o counts_pe.tsv \\
        ${bams.join(' ')}
    """
}

process FEATURE_COUNTS_SE {
    publishDir "${params.outdir}/counts", mode: 'copy', pattern: "counts_se.tsv*"

    input:
    path bams
    path gtf

    output:
    path "counts_se.tsv",         emit: main
    path "counts_se.tsv.summary", emit: summary

    when:
    bams

    script:
    """
    featureCounts -s 2 -T ${task.cpus} \\
        -a ${gtf} \\
        -o counts_se.tsv \\
        ${bams.join(' ')}
    """
}

process MERGE_COUNTS {
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path pe_counts
    path pe_summary
    path se_counts
    path se_summary

    output:
    path "counts.tsv",          emit: main
    path "counts.tsv.summary",  emit: summary

    script:
    // Create dummy files if real ones don't exist to prevent python errors
    if (!pe_counts.exists()) pe_counts.touch()
    if (!pe_summary.exists()) pe_summary.touch()
    if (!se_counts.exists()) se_counts.touch()
    if (!se_summary.exists()) se_summary.touch()
    """
    #!/usr/bin/env python
    import pandas as pd
    import os

    pe_file = '${pe_counts}'
    se_file = '${se_counts}'

    # Check for header and at least one line of data
    pe_has_data = os.path.getsize(pe_file) > 0 and open(pe_file).read().count('\\n') > 1
    se_has_data = os.path.getsize(se_file) > 0 and open(se_file).read().count('\\n') > 1

    # Write a dummy header for featureCounts files
    header = "Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength"
    if pe_has_data:
        df_pe = pd.read_csv(pe_file, sep='\\t', comment='#')
    if se_has_data:
        df_se = pd.read_csv(se_file, sep='\\t', comment='#')

    if pe_has_data and se_has_data:
        merge_cols = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']
        df_final = pd.merge(df_pe, df_se, on=merge_cols)
        df_final.to_csv("counts.tsv", sep='\\t', index=False)
    elif pe_has_data:
        df_pe.to_csv("counts.tsv", sep='\\t', index=False)
    elif se_has_data:
        df_se.to_csv("counts.tsv", sep='\\t', index=False)
    else:
        with open("counts.tsv", 'w') as f:
            f.write(f"{header}\\n")

    # Combine summary files
    with open('counts.tsv.summary', 'w') as outfile:
        if os.path.getsize('${pe_summary}') > 0:
            with open('${pe_summary}') as infile:
                outfile.write(infile.read())
        if os.path.getsize('${se_summary}') > 0:
            with open('${se_summary}') as infile:
                outfile.write("\\n")
                outfile.write(infile.read())
    """
}


process MULTIQC {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path reports

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}
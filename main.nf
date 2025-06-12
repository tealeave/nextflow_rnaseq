#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ========================================================================================
//                          P I P E L I N E  S U M M A R Y
// ========================================================================================
// This pipeline performs RNA-seq analysis with step-wise execution control.
// Use the --step parameter to control how far the pipeline runs:
//   --step qc_raw: Generates FastQC reports on raw reads.
//   --step trim:   Trims reads and generates QC reports on trimmed reads.
//   --step align:  Trims, aligns reads, and stops.
//   --step full:   Runs the complete analysis including feature counting (default).
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
            
            // Set default trim_args if not provided
            if (!meta.trim_args?.trim()) {
                meta.trim_args = params.trim_args_default
            }

            def file_list = [ file(row.fastq_1) ]
            if (row.fastq_2?.trim()) {
                file_list.add(file(row.fastq_2))
            }
            return [meta, file_list]
        }
        .ifEmpty { error "Samplesheet is empty or not found: ${params.samplesheet_file}" }
        .set { ch_raw_reads }

    // --- Step 1: Raw Read QC ---
    if (params.step == 'qc_raw') {
        log.info "Running in QC-only mode (on raw reads)."
        FASTQC_RAW(ch_raw_reads)
        MULTIQC(FASTQC_RAW.out.report.collect().ifEmpty { [] })
        return // Stop workflow here
    }

    // --- Step 2: Read Trimming ---
    ch_raw_reads
        .branch { meta, files ->
            single: meta.trim_type == 'single'
            paired: meta.trim_type == 'paired'
            other:  true
        }
        .set { ch_branched_reads }

    ch_branched_reads.other.count().subscribe { 
        if (it > 0) error "${it} samples found with 'trim_type' other than 'single' or 'paired'." 
    }

    CAT_SE_READS(ch_branched_reads.single)
    TRIM_GALORE_SE(CAT_SE_READS.out.reads)
    TRIM_GALORE_PE(ch_branched_reads.paired)

    ch_trimmed_reads = TRIM_GALORE_SE.out.trimmed_reads.mix(TRIM_GALORE_PE.out.trimmed_reads)
    ch_trim_reports  = TRIM_GALORE_SE.out.report.mix(TRIM_GALORE_PE.out.report)

    FASTQC_TRIM(ch_trimmed_reads)

    if (params.step == 'trim') {
        log.info "Running up to Trimming step."
        def mqc_trim_inputs = Channel.empty().mix(ch_trim_reports).mix(FASTQC_TRIM.out.report).collect()
        MULTIQC(mqc_trim_inputs)
        return // Stop workflow here
    }

    // --- Step 3: Alignment ---
    HISAT2_INDEX(file(params.genome_fasta))
    HISAT2_ALIGN(ch_trimmed_reads, HISAT2_INDEX.out.index_info)

    if (params.step == 'align') {
        log.info "Running up to Alignment step."
        def mqc_align_inputs = Channel.empty()
            .mix(ch_trim_reports)
            .mix(FASTQC_TRIM.out.report)
            .mix(HISAT2_ALIGN.out.summary)
            .collect()
        MULTIQC(mqc_align_inputs)
        return // Stop workflow here
    }

    // --- Step 4: Full Pipeline (Quantification & Final QC) ---
    log.info "Running the full pipeline."

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

    def multiqc_full_inputs = Channel.empty()
        .mix(ch_trim_reports)
        .mix(FASTQC_TRIM.out.report)
        .mix(HISAT2_ALIGN.out.summary)
        .mix(MERGE_COUNTS.out.summary)
        .collect()
    MULTIQC(multiqc_full_inputs)
}

// ========================================================================================
//                        P R O C E S S  D E F I N I T I O N S
// ========================================================================================

// Two separate FASTQC processes for raw and trimmed reads to publish to different dirs
process FASTQC_RAW {
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

process FASTQC_TRIM {
    tag "$meta.id"
    publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy', pattern: '*.{html,zip}'

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
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.fq.gz"), emit: reads

    script:
    """
    cat ${reads.join(' ')} > ${meta.id}.fq.gz
    """
}

process TRIM_GALORE_PE {
    tag "$meta.id"
    // NEW: Organized output and move mode
    publishDir path: "${params.outdir}/trimmed_reads/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_val_*.fq.gz"), emit: trimmed_reads
    path("*_trimming_report.txt"),                  emit: report

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
    // NEW: Organized output and move mode
    publishDir path: "${params.outdir}/trimmed_reads/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fq.gz"), emit: trimmed_reads
    path("*_trimming_report.txt"),                     emit: report

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
    // NEW: publishDir with move mode
    publishDir "${params.outdir}/genome_index", mode: 'copy'

    input:
    path genome_fasta

    output:
    tuple path("hisat2_index"), val(genome_fasta.baseName), emit: index_info

    script:
    """
    mkdir hisat2_index
    hisat2-build -p ${task.cpus} ${genome_fasta} hisat2_index/${genome_fasta.baseName}
    """
}

process HISAT2_ALIGN {
    tag "$meta.id"
    // NEW: Organized output and move mode
    publishDir path: "${params.outdir}/alignment/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    tuple path(index_dir), val(index_prefix)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    path("${meta.id}.summary.txt"),          emit: summary
    path("${meta.id}.bam.bai"),              emit: bai

    script:
    def read_inputs = (meta.trim_type == 'paired') ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads}"
    """
    hisat2 -p ${task.cpus} \\
        --dta \\
        --rna-strandness RF \\
        --summary-file ${meta.id}.summary.txt \\
        -x ${index_dir}/${index_prefix} \\
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
    bams.size() > 0

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
    bams.size() > 0

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
    """
    #!/usr/bin/env python
    import pandas as pd
    import os

    # This Python script remains the same as your original version
    # It correctly handles merging of PE and SE counts
    pe_file = '${pe_counts}'
    se_file = '${se_counts}'

    pe_has_data = os.path.getsize(pe_file) > 0 and open(pe_file).read().count('\\n') > 1
    se_has_data = os.path.getsize(se_file) > 0 and open(se_file).read().count('\\n') > 1

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
    path "multiqc_report.html"

    script:
    """
    multiqc . -f
    """
}

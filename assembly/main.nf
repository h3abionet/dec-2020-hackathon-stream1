#!/usr/bin/env nextflow

// Get sample info from sample sheet
// Minimum information that are needed in the sample sheet are SampleID, forward and reverse read file location
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def fastq_r1_file = file(row['FastqR1'])
            def fastq_r2_file = file(row['FastqR2'])
            return [ sample_id, fastq_r1_file, fastq_r2_file ]
        }.into{samples_1; samples_2}


process print_sample_info {
    tag { "${sample_id}" }
    echo true
    input:
    set val(sample_id), file(fastq_r1_file), file(fastq_r2_file) from samples_1
    script:
    """
    printf "[sample_info] sample: ${sample_id}\tFastqR1: ${fastq_r1_file}\tFastqR2: ${fastq_r2_file}\n"
    """
}

process log_version_sga {
    tag { "${params.project_name}.logVersionSga" }
    echo true
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'bwa_samtools'

    output:
    file("tool.sga.version") into tool_version_sga

    script:
    """
    sga --version > tool.sga.version
    """
}

process run_sga_preprocess {
    tag { "${params.project_name}.${sample_id}.rSgaPreProcess" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file(fastq_r1_file), file(fastq_r2_file) from samples_2

    output:
    set val(sample_id), file("${sample_id}.fastq") into pp_fastq

    script:
    """
    sga preprocess \
    -o ${sample_id}.fastq \
    --pe-mode 1 \
    ${fastq_r1_file} \
    ${fastq_r2_file}
    """
}

process run_sga_index {
    tag { "${params.project_name}.${sample_id}.rSgaIndex" }
    memory { 64.GB * task.attempt }
    cpus { "${params.sga_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file(fastq_file) from fastq

    output:
    set val(sample_id), file("${sample_id}.*")  into index

    script:
    """
    sga index \
    -a ropebwt \
    -no-reverse \
    -t ${params.sga_threads} \
    ${fastq_file}
    """
}

process run_sga_correct {
    tag { "${params.project_name}.${sample_id}.rSgaCorrect" }
    memory { 64.GB * task.attempt }
    cpus { "${params.sga_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'symlink', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file("${sample_id}.fastq") from index

    output:
    set val(sample_id), file("${sample_id}.*")  into correct

    script:
    """
    sga correct \
    -k ${k} \
    --learn \
    -t ${params.sga_threads} \
    -o ${sample_id}.correct.fastq \
    ${sample_id}.fastq
    """
}

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

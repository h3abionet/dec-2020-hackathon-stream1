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
        }.set{samples}

process sga_log_version {
    tag { "${params.project_name}.sgaLogVersion" }
    echo true
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'bwa_samtools'

    output:
    file("sga.version") into log_version

    script:
    """
    sga --version > sga.version
    """
}

process sga_preprocess {
    tag { "${params.project_name}.${sample_id}.sgaPreProcess" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file(fastq_r1_file), file(fastq_r2_file) from samples

    output:
    set val(sample_id), file("${sample_id}.fastq") into preprocess

    script:
    """
    sga preprocess \
    -o ${sample_id}.fastq \
    --pe-mode 1 \
    ${fastq_r1_file} \
    ${fastq_r2_file}
    """
}

process sga_index {
    tag { "${params.project_name}.${sample_id}.sgaIndex" }
    cpus { "${params.sga_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file(fastq_file) from preprocess

    output:
    set val(sample_id), file(fastq_file), file("${sample_id}.bwt"), file("${sample_id}.sai") into index

    script:
    """
    sga index \
    -a ropebwt \
    --no-reverse \
    -t ${params.sga_threads} \
    ${fastq_file}
    """
}

process sga_correct {
    tag { "${params.project_name}.${sample_id}.sgaCorrect" }
    cpus { "${params.sga_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file("${sample_id}.fastq"), file("${sample_id}.bwt"), file("${sample_id}.sai") from index

    output:
    set val(sample_id), file("${sample_id}.correct.fastq")  into correct

    script:
    """
    sga correct \
    -k ${params.kmer_size} \
    --learn \
    -t ${params.sga_threads} \
    -o ${sample_id}.correct.fastq \
    ${sample_id}.fastq
    """
}

process sga_index_on_correct {
    tag { "${params.project_name}.${sample_id}.sgaIndexOnCorrect" }
    cpus { "${params.sga_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file("${sample_id}.correct.fastq") from correct

    output:
    set val(sample_id), file("${sample_id}.correct.fastq"), file("${sample_id}.correct.bwt"), file("${sample_id}.correct.sai"), file("${sample_id}.correct.rbwt"), file("${sample_id}.correct.rsai") into index_on_correct

    script:
    """
    sga index \
    -a ropebwt \
    -t ${params.sga_threads} \
    ${sample_id}.correct.fastq
    """
}

process sga_filter {
    tag { "${params.project_name}.${sample_id}.sgaFilter" }
    cpus { "${params.sga_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file("${sample_id}.correct.fastq"), file("${sample_id}.correct.bwt"), file("${sample_id}.correct.sai"), file("${sample_id}.correct.rbwt"), file("${sample_id}.correct.rsai") from index_on_correct

    output:
    set val(sample_id), file("${sample_id}.correct.filter.pass.bwt"), file("${sample_id}.correct.filter.pass.sai"), file("${sample_id}.correct.filter.pass.rbwt"), file("${sample_id}.correct.filter.pass.rsai"), file("${sample_id}.correct.filter.pass.fa") into filter

    script:
    """
    sga filter \
    -x 2 \
    -t ${params.sga_threads} \
    -o ${sample_id}.correct.filter.pass.fa \
    ${sample_id}.correct.fastq
    """
}

process sga_fm_merge {
    tag { "${params.project_name}.${sample_id}.sgaFmMerge" }
    cpus { "${params.sga_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file("${sample_id}.correct.filter.pass.bwt"), file("${sample_id}.correct.filter.pass.sai"), file("${sample_id}.correct.filter.pass.rbwt"), file("${sample_id}.correct.filter.pass.rsai"), file("${sample_id}.correct.filter.pass.fa") from filter

    output:
    set val(sample_id), file("${sample_id}.correct.filter.pass.merged.fa") into fm_merge

    script:
    """
    sga fm-merge \
    -m ${params.min_overlap} \
    -t ${params.sga_threads} \
    -o ${sample_id}.correct.filter.pass.merged.fa \
    ${sample_id}.correct.filter.pass.fa
    """
}

process sga_index_on_fm_merge {
    tag { "${params.project_name}.${sample_id}.sgaIndexOnFmMerge" }
    cpus { "${params.sga_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file("${sample_id}.correct.filter.pass.merged.fa") from fm_merge

    output:
    set val(sample_id), file("${sample_id}.correct.filter.pass.merged.fa"), file("${sample_id}.correct.filter.pass.merged.bwt"), file("${sample_id}.correct.filter.pass.merged.sai"),  file("${sample_id}.correct.filter.pass.merged.rbwt"), file("${sample_id}.correct.filter.pass.merged.rsai") into index_on_fm_merge

    script:
    """
    sga index \
    -d 20000000 \
    -t ${params.sga_threads} \
    ${sample_id}.correct.filter.pass.merged.fa
    """
}

process sga_rmdup {
    tag { "${params.project_name}.${sample_id}.sgaRmdup" }
    cpus { "${params.sga_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file("${sample_id}.correct.filter.pass.merged.fa"), file("${sample_id}.correct.filter.pass.merged.bwt"), file("${sample_id}.correct.filter.pass.merged.sai"), file("${sample_id}.correct.filter.pass.merged.rbwt"), file("${sample_id}.correct.filter.pass.merged.rsai") from index_on_fm_merge 

    output:
    set val(sample_id), file("${sample_id}.correct.filter.pass.merged.rmdup.fa"), file("${sample_id}.correct.filter.pass.merged.rmdup.bwt"), file("${sample_id}.correct.filter.pass.merged.rmdup.sai"),  file("${sample_id}.correct.filter.pass.merged.rmdup.rbwt"), file("${sample_id}.correct.filter.pass.merged.rmdup.rsai") into rmdup

    script:
    """
    sga rmdup \
    -t ${params.sga_threads} \
    -o ${sample_id}.correct.filter.pass.merged.rmdup.fa \
    ${sample_id}.correct.filter.pass.merged.fa
    """
}

process sga_overlap {
    tag { "${params.project_name}.${sample_id}.sgaOverlap" }
    cpus { "${params.sga_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file("${sample_id}.correct.filter.pass.merged.rmdup.fa"), file("${sample_id}.correct.filter.pass.merged.rmdup.bwt"), file("${sample_id}.correct.filter.pass.merged.rmdup.sai"),  file("${sample_id}.correct.filter.pass.merged.rmdup.rbwt"), file("${sample_id}.correct.filter.pass.merged.rmdup.rsai") from rmdup
  
    output:
    set val(sample_id), file("${sample_id}.correct.filter.pass.merged.rmdup.asqg.gz") into overlap

    script:
    """
    sga overlap \
    -m ${params.min_overlap} \
    -t ${params.sga_threads} \
    ${sample_id}.correct.filter.pass.merged.rmdup.fa
    """
}

process sga_assemble {
    tag { "${params.project_name}.${sample_id}.sgaAssemble" }
    cpus { "${params.sga_threads}" }
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy', overwrite: false
    label 'sga'

    input:
    set val(sample_id), file("${sample_id}.correct.filter.pass.merged.rmdup.asqg.gz") from overlap

    output:
    set val(sample_id), file("${sample_id}.assemble-contigs.fa"), file("${sample_id}.assemble-graph.asqg.gz"), file("${sample_id}.assemble-variants.fa") into assemble

    script:
    """
    sga assemble \
    -m ${params.assembly_min_overlap} \
    -l ${params.min_branch_length} \
    -o ${sample_id}.assemble \
    ${sample_id}.correct.filter.pass.merged.rmdup.asqg.gz
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

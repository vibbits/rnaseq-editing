// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)
    
    if (params.enable_aks) {
       pod nodeSelector: 'agentpool=cpumem'
    }
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/trimmomatic:0.39--0"
    } else {
        container "quay.io/biocontainers/trimmomatic:0.39--0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq.gz")    , emit: reads
    //tuple val(meta), path("*report.txt"), emit: log
    path "*.version.txt"                , emit: version

    //tuple val(meta), path("*.html"), emit: html optional true
    //tuple val(meta), path("*.zip") , emit: zip optional true

    script:
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (meta.single_end) cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 4) cores = 6 
    }

    // Clipping presets have to be evaluated in the context of SE/PE
    def c_lead   = params.leading > 0             ? "LEADING:${params.leading}"                         : ''
    def c_trail  = params.trailing > 0             ? "TRAILING:${params.trailing}"                         : ''
    def c_extra  = params.trimextra ? "${params.trimextra}" : ''
    def c_adapter = params.adapter_file_path + params.adapter_file 
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        [ ! -f $c_adapter ] && ln -s ${params.adapter_file} $c_adapter 
        trimmomatic SE\\
            $options.args \\
            -threads $cores \\
            ${prefix}_1.fastq.gz \\
            ${prefix}_1.fq.gz \\
            ILLUMINACLIP:${params.adapter_file}:2:30:10:1:true \\
            $c_lead \\
            $c_trail \\
            $c_extra

        echo '0.39' > ${software}.version.txt
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        [ ! -f $c_adapter ] && ln -s ${params.adapter_file} $c_adapter 
        trimmomatic PE\\
            $options.args \\
            -threads $cores \\
            ${prefix}_1.fastq.gz \\
            ${prefix}_2.fastq.gz \\
            ${prefix}_1.fq.gz \\
            ${prefix}_U_1.fastq.gz \\
            ${prefix}_2.fq.gz \\
            ${prefix}_U_2.fastq.gz \\
            ILLUMINACLIP:${params.adapter_file}:2:30:10:1:true \\
            $c_lead \\
            $c_trail \\
            $c_extra
        echo '0.39' > ${software}.version.txt
        """
    }
}

// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process GET_CHROM_SIZES {
    label 'process_low'
    tag "$fasta"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'genome', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::samtools=1.10" : null)
    if (params.enable_aks) {
       pod nodeSelector: 'agentpool=cpumem'
    }   
 
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"
    } else {
        container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    }

    input:
    path fasta

    output:
    path '*.sizes'      , emit: sizes
    path '*.fai'        , emit: fai
    path "*.version.txt", emit: version

    script:
    def software = 'samtools'
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}

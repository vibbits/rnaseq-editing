// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RDDPRED_PREPAREGROUP {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    
    file(conditions)  
    val(group)

    output:
    path '*.txt', emit: input_conditions 

    script:
    """
    mkdir -p tmp/conds
    ## read samplesheet and extract groups defined by params   
    ## check sample names of files 
    grep $group $conditions > tmp/conds/$groups.txt
    ## move files with group to output
    mv tmp/conds/*.txt . 
    """
}

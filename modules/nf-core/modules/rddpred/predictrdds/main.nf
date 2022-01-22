// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RDDPRED_PREDICTRDDS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda (params.enable_conda ? 'bioconda::star=2.6.1d' : null)
    if (params.enable_aks) {
       pod nodeSelector: 'agentpool=cpumem'
    }


    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'rddpred-github:1.1.3'
    } else {
        container 'docker.io/vibbioinfocore/rddpred:1.1.4'
    }

    input:
    tuple val(groups), path(groupfile)
    file params.fasta 

    output:
    path '*.version.txt', emit: version
    path '*.Log', emit: logs
    path '*.Step2*.txt', emit: step2
    path '*.Step5.Results.Summary.txt', emit: summary
    path '*.ModelDir', emit: modeldir
    path '*.Prediction.ResultList.txt', emit: resultlist
    path '*.RDDpred.results_report.txt', emit: resultreport
    path '*.RDD.RawList.txt', emit: rawlist
    path '*.VafList.txt', emit: vaflist
    
    script:
    // Calculate number of --cores for RDDpred based on value of task.cpus
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = task.cpus 
    //if (task.cpus) {
    //    cores = (task.cpus as int) - 4
    //    if (cores < 1) cores = 1
    //    if (cores > 4) cores = 20 
    //}

    def software   = getSoftwareName(task.process)
    def positive = params.rddpred_databases_file_path + params.pos_site_list
    def negative = params.rddpred_databases_file_path + params.neg_site_list
    """
    [ ! -f $negative ] && ln -s $negative ${params.neg_site_list}
    [ ! -f $positive ] && ln -s $positive ${params.pos_site_list}
    python /code/RDDpred.py \\
        -rsf ${params.fasta} \\
        -rbl $groupfile  \\
        -pni $cores \\
        -ops $groups \\
        -psl $positive \\
        -nsl $negative \\
        $options.args

    echo '1.1.3' > ${software}.version.txt
    rm $groupfile
    """
}

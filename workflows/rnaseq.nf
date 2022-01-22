/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def valid_params = [
    aligners       : ['star_salmon', 'star_rsem', 'hisat2'],
    pseudoaligners : ['salmon'],
    rseqc_modules  : ['bam_stat', 'inner_distance', 'infer_experiment', 'junction_annotation', 'junction_saturation', 'read_distribution', 'read_duplication']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnaseq.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta, params.transcript_fasta, params.additional_fasta,
    params.gtf, params.gff, params.gene_bed,
    params.ribo_database_manifest, params.splicesites,
    params.star_index, params.hisat2_index, params.rsem_index, params.salmon_index
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check rRNA databases for sortmerna
ch_ribo_db = file(params.ribo_database_manifest)
if (ch_ribo_db.isEmpty()) {exit 1, "File ${ch_ribo_db.getName()} is empty!"}

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_alignment) { prepareToolIndices << params.aligner        }
//if (params.pseudo_aligner)  { prepareToolIndices << params.pseudo_aligner }

// Get RSeqC modules to run
def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// Header files for MultiQC
ch_pca_header_multiqc        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_clustering_header_multiqc = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
ch_biotypes_header_multiqc   = file("$projectDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]

def multiqc_options         = modules['multiqc']
multiqc_options.args       += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''
if (params.skip_alignment)  { multiqc_options['publish_dir'] = '' }

def deseq2_qc_options                 = modules['deseq2_qc']
deseq2_qc_options.args               += params.deseq2_vst ? Utils.joinModuleArgs(['--vst TRUE']) : ''
def deseq2_qc_salmon_options          = deseq2_qc_options.clone()
deseq2_qc_salmon_options.publish_dir  = "salmon/deseq2_qc"

include { BEDTOOLS_GENOMECOV                 } from '../modules/local/bedtools_genomecov'          addParams( options: modules['bedtools_genomecov']                     )
include { DESEQ2_QC as DESEQ2_QC_STAR_SALMON } from '../modules/local/deseq2_qc'                   addParams( options: deseq2_qc_options, multiqc_label: 'star_salmon'   )
include { DESEQ2_QC as DESEQ2_QC_RSEM        } from '../modules/local/deseq2_qc'                   addParams( options: deseq2_qc_options, multiqc_label: 'star_rsem'     )
include { DESEQ2_QC as DESEQ2_QC_SALMON      } from '../modules/local/deseq2_qc'                   addParams( options: deseq2_qc_salmon_options, multiqc_label: 'salmon' )
include { DUPRADAR                           } from '../modules/local/dupradar'                    addParams( options: modules['dupradar']                               )
include { GET_SOFTWARE_VERSIONS              } from '../modules/local/get_software_versions'       addParams( options: [publish_files : ['tsv':'']]                      )
include { MULTIQC                            } from '../modules/local/multiqc'                     addParams( options: multiqc_options                                   )
include { MULTIQC_CUSTOM_BIOTYPE             } from '../modules/local/multiqc_custom_biotype'      addParams( options: modules['multiqc_custom_biotype']                 )
include { MULTIQC_CUSTOM_FAIL_MAPPED         } from '../modules/local/multiqc_custom_fail_mapped'  addParams( options: [publish_files: false]                            )
include { MULTIQC_CUSTOM_STRAND_CHECK        } from '../modules/local/multiqc_custom_strand_check' addParams( options: [publish_files: false]                            )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
def gffread_options         = modules['gffread']
if (!params.save_reference) { gffread_options['publish_files'] = false }

def star_genomegenerate_options = modules['star_genomegenerate']
if (!params.save_reference)     { star_genomegenerate_options['publish_files'] = false }

def hisat2_build_options    = modules['hisat2_build']
if (!params.save_reference) { hisat2_build_options['publish_files'] = false }

def rsem_preparereference_options = modules['rsem_preparereference']
if (!params.save_reference)       { rsem_preparereference_options['publish_files'] = false }

def rsem_calculateexpression_options = modules['rsem_calculateexpression']
if (params.save_align_intermeds)     { rsem_calculateexpression_options.publish_files.put('bam','') }

def salmon_index_options     = modules['salmon_index']
salmon_index_options.args   += params.gencode ? Utils.joinModuleArgs(['--gencode']) : ''
if (!params.save_reference)  { salmon_index_options['publish_files'] = false }

def salmon_quant_options   = modules['salmon_quant']
salmon_quant_options.args += params.save_unaligned ? Utils.joinModuleArgs(['--writeUnmappedNames']) : ''

def samtools_sort_genome_options  = modules['samtools_sort_genome']
def samtools_index_genome_options = modules['samtools_index_genome']
samtools_index_genome_options.args += params.bam_csi_index ? Utils.joinModuleArgs(['-c']) : ''
if (['star_salmon','hisat2'].contains(params.aligner)) {
    if (params.save_align_intermeds || (!params.with_umi && params.skip_markduplicates)) {
        samtools_sort_genome_options.publish_files.put('bam','')
        samtools_index_genome_options.publish_files.put('bai','')
        samtools_index_genome_options.publish_files.put('csi','')
    }
} else {
    if (params.save_align_intermeds || params.skip_markduplicates) {
        samtools_sort_genome_options.publish_files.put('bam','')
        samtools_index_genome_options.publish_files.put('bai','')
        samtools_index_genome_options.publish_files.put('csi','')
    }
}

include { INPUT_CHECK    } from '../subworkflows/local/input_check'    addParams( options: [:] )
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome' addParams( genome_options: publish_genome_options, index_options: publish_index_options, gffread_options: gffread_options,  star_index_options: star_genomegenerate_options,  hisat2_index_options: hisat2_build_options, rsem_index_options: rsem_preparereference_options, salmon_index_options: salmon_index_options )
include { QUANTIFY_RSEM  } from '../subworkflows/local/quantify_rsem'  addParams( calculateexpression_options: rsem_calculateexpression_options, samtools_sort_options: samtools_sort_genome_options, samtools_index_options: samtools_index_genome_options, samtools_stats_options: samtools_index_genome_options, merge_counts_options: modules['rsem_merge_counts'] )
include { QUANTIFY_SALMON as QUANTIFY_STAR_SALMON } from '../subworkflows/local/quantify_salmon'    addParams( genome_options: publish_genome_options, tximport_options: modules['star_salmon_tximport'], salmon_quant_options: modules['star_salmon_quant'], merge_counts_options: modules['star_salmon_merge_counts'] )
include { QUANTIFY_SALMON as QUANTIFY_SALMON      } from '../subworkflows/local/quantify_salmon'    addParams( genome_options: publish_genome_options, tximport_options: modules['salmon_tximport'], salmon_quant_options: salmon_quant_options, merge_counts_options: modules['salmon_merge_counts'] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
def fastqc_options = modules['fastqc']
def fastqc_trim_options = modules['fastqc'].clone()
fastqc_trim_options.publish_dir = "fastqc_trim"
fastqc_trim_options.suffix = 'trimmed'

def rddpred_options          = modules['rddpred']
rddpred_options.publish_dir = "rddpred"

def cat_fastq_options          = modules['cat_fastq']
if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

def sortmerna_options           = modules['sortmerna']
if (params.save_non_ribo_reads) { sortmerna_options.publish_files.put('fastq.gz','') }

def stringtie_options   = modules['stringtie']
stringtie_options.args += params.stringtie_ignore_gtf ? '' : Utils.joinModuleArgs(['-e'])

def subread_featurecounts_options  = modules['subread_featurecounts']
def biotype                        = params.gencode ? "gene_type" : params.featurecounts_group_type
subread_featurecounts_options.args += Utils.joinModuleArgs(["-g $biotype", "-t $params.featurecounts_feature_type"])

include { RDDPRED_PREDICTRDDS   } from '../modules/nf-core/modules/rddpred/predictrdds/main'   addParams( options: rddpred_options                            )
include { CAT_FASTQ             } from '../modules/nf-core/modules/cat/fastq/main'             addParams( options: cat_fastq_options                            )
include { SAMTOOLS_SORT         } from '../modules/nf-core/modules/samtools/sort/main'         addParams( options: modules['umitools_dedup_transcriptome_sort'] )
include { PRESEQ_LCEXTRAP       } from '../modules/nf-core/modules/preseq/lcextrap/main'       addParams( options: modules['preseq_lcextrap']                   )
include { QUALIMAP_RNASEQ       } from '../modules/nf-core/modules/qualimap/rnaseq/main'       addParams( options: modules['qualimap_rnaseq']                   )
include { SORTMERNA             } from '../modules/nf-core/modules/sortmerna/main'             addParams( options: sortmerna_options                            )
include { STRINGTIE             } from '../modules/nf-core/modules/stringtie/stringtie/main'   addParams( options: stringtie_options                            )
include { SUBREAD_FEATURECOUNTS } from '../modules/nf-core/modules/subread/featurecounts/main' addParams( options: subread_featurecounts_options                )
include { FASTQC as FASTQC_TRIM } from '../modules/nf-core/modules/fastqc/main'                addParams( options: fastqc_trim_options     )

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
def umitools_extract_options    = modules['umitools_extract']

umitools_extract_options.args  += params.umitools_extract_method ? Utils.joinModuleArgs(["--extract-method=${params.umitools_extract_method}"]) : ''
umitools_extract_options.args  += params.umitools_bc_pattern     ? Utils.joinModuleArgs(["--bc-pattern='${params.umitools_bc_pattern}'"])       : ''
if (params.save_umi_intermeds)  { umitools_extract_options.publish_files.put('fastq.gz','') }

def trimmomatic_options = modules['trimmomatic']
if (params.save_trimmed) {trimmomatic_options.publish_files.put('fq.gz','') }

def trimgalore_options    = modules['trimgalore']
trimgalore_options.args  += params.trim_nextseq > 0 ? Utils.joinModuleArgs(["--nextseq ${params.trim_nextseq}"]) : ''
if (params.save_trimmed)  { trimgalore_options.publish_files.put('fq.gz','') }

def star_align_options            = modules['star_align']
star_align_options.args          += params.save_unaligned ? Utils.joinModuleArgs(['--outReadsUnmapped Fastx']) : ''
if (params.save_align_intermeds)  { star_align_options.publish_files.put('bam','') }
if (params.save_unaligned)        { star_align_options.publish_files.put('fastq.gz','unmapped') }

def hisat2_align_options         = modules['hisat2_align']
if (params.save_align_intermeds) { hisat2_align_options.publish_files.put('bam','') }
if (params.save_unaligned)       { hisat2_align_options.publish_files.put('fastq.gz','unmapped') }

def picard_markduplicates_samtools   = modules['picard_markduplicates_samtools']
picard_markduplicates_samtools.args += params.bam_csi_index ? Utils.joinModuleArgs(['-c']) : ''

def umitools_dedup_genome_options            = modules['umitools_dedup_genome']
def umitools_dedup_genome_samtools_options   = modules['umitools_dedup_genome_samtools']
umitools_dedup_genome_samtools_options.args += params.bam_csi_index ? Utils.joinModuleArgs(['-c']) : ''
if (['star_salmon','hisat2'].contains(params.aligner)) {
    if (params.save_align_intermeds || params.skip_markduplicates || params.save_umi_intermeds) {
        umitools_dedup_genome_options.publish_files.put('bam','')
        umitools_dedup_genome_samtools_options.publish_files.put('bai','')
        umitools_dedup_genome_samtools_options.publish_files.put('csi','')
    }
}

include { FASTQC_UMITOOLS_TRIMMOMATIC } from '../subworkflows/nf-core/fastqc_umitools_trimmomatic' addParams( fastqc_options: fastqc_options, umitools_options: umitools_extract_options, trimmomatic_options: trimmomatic_options ) 

include { FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastqc_umitools_trimgalore' addParams( fastqc_options: fastqc_options, umitools_options: umitools_extract_options, trimgalore_options: trimgalore_options )

include { ALIGN_STAR                 } from '../subworkflows/nf-core/align_star'                 addParams( align_options: star_align_options, samtools_sort_options: samtools_sort_genome_options, samtools_index_options: samtools_index_genome_options, samtools_stats_options: samtools_index_genome_options   )
include { ALIGN_HISAT2               } from '../subworkflows/nf-core/align_hisat2'               addParams( align_options: hisat2_align_options, samtools_sort_options: samtools_sort_genome_options, samtools_index_options: samtools_index_genome_options, samtools_stats_options: samtools_index_genome_options )
include { BAM_SORT_SAMTOOLS          } from '../subworkflows/nf-core/bam_sort_samtools'          addParams( sort_options: modules['samtools_sort_transcriptome'], index_options: modules['samtools_index_transcriptome'], stats_options: modules['samtools_index_transcriptome']      )
include { MARK_DUPLICATES_PICARD     } from '../subworkflows/nf-core/mark_duplicates_picard'     addParams( markduplicates_options: modules['picard_markduplicates'], samtools_index_options: picard_markduplicates_samtools, samtools_stats_options:  picard_markduplicates_samtools )
include { RSEQC                      } from '../subworkflows/nf-core/rseqc'                      addParams( bamstat_options: modules['rseqc_bamstat'], innerdistance_options: modules['rseqc_innerdistance'], inferexperiment_options: modules['rseqc_inferexperiment'], junctionannotation_options: modules['rseqc_junctionannotation'], junctionsaturation_options: modules['rseqc_junctionsaturation'], readdistribution_options: modules['rseqc_readdistribution'], readduplication_options: modules['rseqc_readduplication'] )
include { DEDUP_UMI_UMITOOLS as DEDUP_UMI_UMITOOLS_GENOME        } from '../subworkflows/nf-core/dedup_umi_umitools' addParams( umitools_options: umitools_dedup_genome_options, samtools_index_options: umitools_dedup_genome_samtools_options, samtools_stats_options: umitools_dedup_genome_samtools_options             )
include { DEDUP_UMI_UMITOOLS as DEDUP_UMI_UMITOOLS_TRANSCRIPTOME } from '../subworkflows/nf-core/dedup_umi_umitools' addParams( umitools_options: modules['umitools_dedup_transcriptome'], samtools_index_options: modules['umitools_dedup_transcriptome'], samtools_stats_options: modules['umitools_dedup_transcriptome'] )
include { BEDGRAPH_TO_BIGWIG as BEDGRAPH_TO_BIGWIG_FORWARD       } from '../subworkflows/nf-core/bedgraph_to_bigwig' addParams( bedclip_options: modules['ucsc_bedclip_forward'], bedgraphtobigwig_options: modules['ucsc_bedgraphtobigwig_forward'] )
include { BEDGRAPH_TO_BIGWIG as BEDGRAPH_TO_BIGWIG_REVERSE       } from '../subworkflows/nf-core/bedgraph_to_bigwig' addParams( bedclip_options: modules['ucsc_bedclip_reverse'], bedgraphtobigwig_options: modules['ucsc_bedgraphtobigwig_reverse'] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report      = []
def pass_percent_mapped = [:]
def fail_percent_mapped = [:]

workflow RNASEQ {

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (
        prepareToolIndices,
        biotype
    )
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.gffread_version.ifEmpty(null))

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_fastq

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    
    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters
    //
    FASTQC_UMITOOLS_TRIMMOMATIC (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.with_umi,
        params.skip_trimming
    )
   ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMMOMATIC.out.fastqc_version.first().ifEmpty(null))
   ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMMOMATIC.out.umitools_version.first().ifEmpty(null))
   ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMMOMATIC.out.trimmomatic_version.first().ifEmpty(null))
   
   ch_trimmed_reads_fastqc = FASTQC_UMITOOLS_TRIMMOMATIC.out.reads
   // input subworkflow alignment with STAR
   ch_trimmed_reads = ch_trimmed_reads_fastqc 
   
   if(!params.skip_trimming) { 
   FASTQC_TRIM (
      ch_trimmed_reads_fastqc
   ) 
   ch_trimmed_fastqc_html_reports = FASTQC_TRIM.out.html
   ch_trimmed_fastqc_zip = FASTQC_TRIM.out.zip
   ch_trimmed_software_versions = FASTQC_TRIM.out.version   
   }

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //
    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_salmon') {
        ALIGN_STAR (
            ch_trimmed_reads,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.gtf
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final
        if (params.bam_csi_index) {
            ch_genome_bam_index  = ALIGN_STAR.out.csi
        }
        ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.star_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.samtools_version.first().ifEmpty(null))
 
    }

    // create bamlists files for RDDPred based on the sample name <sample>_<condition>_<replicate>
    if (!params.skip_rddpred) {
       ch_genome_bam_for_rddpred = Channel.empty()
       ch_genome_bam_for_rddpred = ch_genome_bam.collect().flatten().collate(2)
 
       ch_genome_bam_for_rddpred
          .map((sample_info,path) -> [sample_info.id.split('_')[1], path])
          .groupTuple()
          .collectFile( storeDir: "$projectDir/assets/" )  {
             it -> ["${it[0]}.conditions.txt", it[1].join('\n')]     }

       ch_rddpred_bam_input = Channel.empty()     
       ch_rddpred_bam_input = Channel.fromPath("$projectDir/assets/*.conditions.txt")
         .map { file -> tuple(file.baseName.split('.conditions.txt')[0], file) }
         .view()

       // launch RDDpred

       RDDPRED_PREDICTRDDS (
          ch_rddpred_bam_input,
          PREPARE_GENOME.out.star_index
       )  
    }

    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with RSEM
    //
    ch_rsem_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_rsem') {
        QUANTIFY_RSEM (
            ch_trimmed_reads,
            PREPARE_GENOME.out.rsem_index
        )
        ch_genome_bam        = QUANTIFY_RSEM.out.bam
        ch_genome_bam_index  = QUANTIFY_RSEM.out.bai
        ch_samtools_stats    = QUANTIFY_RSEM.out.stats
        ch_samtools_flagstat = QUANTIFY_RSEM.out.flagstat
        ch_samtools_idxstats = QUANTIFY_RSEM.out.idxstats
        ch_star_multiqc      = QUANTIFY_RSEM.out.logs
        ch_rsem_multiqc      = QUANTIFY_RSEM.out.stat
        if (params.bam_csi_index) {
            ch_genome_bam_index  = QUANTIFY_RSEM.out.csi
        }
        ch_software_versions = ch_software_versions.mix(QUANTIFY_RSEM.out.rsem_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(QUANTIFY_RSEM.out.samtools_version.first().ifEmpty(null))

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_RSEM (
                QUANTIFY_RSEM.out.merged_counts_gene,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_aligner_pca_multiqc        = DESEQ2_QC_RSEM.out.pca_multiqc
            ch_aligner_clustering_multiqc = DESEQ2_QC_RSEM.out.dists_multiqc
            ch_software_versions          = ch_software_versions.mix(DESEQ2_QC_RSEM.out.version.ifEmpty(null))
        }
    }

    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    ch_hisat2_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        ALIGN_HISAT2 (
            ch_trimmed_reads,
            PREPARE_GENOME.out.hisat2_index,
            PREPARE_GENOME.out.splicesites
        )
        ch_genome_bam        = ALIGN_HISAT2.out.bam
        ch_genome_bam_index  = ALIGN_HISAT2.out.bai
        ch_samtools_stats    = ALIGN_HISAT2.out.stats
        ch_samtools_flagstat = ALIGN_HISAT2.out.flagstat
        ch_samtools_idxstats = ALIGN_HISAT2.out.idxstats
        ch_hisat2_multiqc    = ALIGN_HISAT2.out.summary
        if (params.bam_csi_index) {
            ch_genome_bam_index  = ALIGN_HISAT2.out.csi
        }
        ch_software_versions = ch_software_versions.mix(ALIGN_HISAT2.out.hisat2_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_HISAT2.out.samtools_version.first().ifEmpty(null))

        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //
        if (params.with_umi) {
            DEDUP_UMI_UMITOOLS_GENOME (
                ch_genome_bam.join(ch_genome_bam_index, by: [0])
            )
            ch_genome_bam        = DEDUP_UMI_UMITOOLS_GENOME.out.bam
            ch_genome_bam_index  = DEDUP_UMI_UMITOOLS_GENOME.out.bai
            ch_samtools_stats    = DEDUP_UMI_UMITOOLS_GENOME.out.stats
            ch_samtools_flagstat = DEDUP_UMI_UMITOOLS_GENOME.out.flagstat
            ch_samtools_idxstats = DEDUP_UMI_UMITOOLS_GENOME.out.idxstats
            if (params.bam_csi_index) {
                ch_genome_bam_index  = DEDUP_UMI_UMITOOLS_GENOME.out.csi
            }
        }
    }

    //
    // Filter channels to get samples that passed STAR minimum mapping percentage
    //
    ch_fail_mapping_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner.contains('star')) {
        ch_star_multiqc
            .map { meta, align_log -> [ meta ] + WorkflowRnaseq.getStarPercentMapped(params, align_log) }
            .set { ch_percent_mapped }

        ch_genome_bam
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam }

        ch_genome_bam_index
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam_index }

        ch_percent_mapped
            .branch { meta, mapped, pass ->
                pass: pass
                    pass_percent_mapped[meta.id] = mapped
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    fail_percent_mapped[meta.id] = mapped
                    return [ "$meta.id\t$mapped" ]
            }
            .set { ch_pass_fail_mapped }

        MULTIQC_CUSTOM_FAIL_MAPPED (
            ch_pass_fail_mapped.fail.collect()
        )
        .set { ch_fail_mapping_multiqc }
    }

    //
    // MODULE: Downstream QC steps
    //
    ch_qualimap_multiqc           = Channel.empty()
    ch_dupradar_multiqc           = Channel.empty()
    ch_bamstat_multiqc            = Channel.empty()
    ch_inferexperiment_multiqc    = Channel.empty()
    ch_innerdistance_multiqc      = Channel.empty()
    ch_junctionannotation_multiqc = Channel.empty()
    ch_junctionsaturation_multiqc = Channel.empty()
    ch_readdistribution_multiqc   = Channel.empty()
    ch_readduplication_multiqc    = Channel.empty()
    ch_fail_strand_multiqc        = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc) {
        if (!params.skip_qualimap) {
            QUALIMAP_RNASEQ (
                ch_genome_bam,
                PREPARE_GENOME.out.gtf
            )
            ch_qualimap_multiqc  = QUALIMAP_RNASEQ.out.results
            ch_software_versions = ch_software_versions.mix(QUALIMAP_RNASEQ.out.version.first().ifEmpty(null))
        }
        if (!params.skip_dupradar) {
            DUPRADAR (
                ch_genome_bam,
                PREPARE_GENOME.out.gtf
            )
            ch_dupradar_multiqc  = DUPRADAR.out.multiqc
            ch_software_versions = ch_software_versions.mix(DUPRADAR.out.version.first().ifEmpty(null))
        }
        if (!params.skip_rseqc && rseqc_modules.size() > 0) {
            RSEQC (
                ch_genome_bam,
                PREPARE_GENOME.out.gene_bed,
                rseqc_modules
            )
            ch_bamstat_multiqc            = RSEQC.out.bamstat_txt
            ch_inferexperiment_multiqc    = RSEQC.out.inferexperiment_txt
            ch_innerdistance_multiqc      = RSEQC.out.innerdistance_freq
            ch_junctionannotation_multiqc = RSEQC.out.junctionannotation_log
            ch_junctionsaturation_multiqc = RSEQC.out.junctionsaturation_rscript
            ch_readdistribution_multiqc   = RSEQC.out.readdistribution_txt
            ch_readduplication_multiqc    = RSEQC.out.readduplication_pos_xls
            ch_software_versions          = ch_software_versions.mix(RSEQC.out.version.first().ifEmpty(null))

            ch_inferexperiment_multiqc
                .map { meta, strand_log -> [ meta ] + WorkflowRnaseq.getInferexperimentStrandedness(strand_log, 30) }
                .filter { it[0].strandedness != it[1] }
                .map { meta, strandedness, sense, antisense, undetermined ->
                    [ "$meta.id\t$meta.strandedness\t$strandedness\t$sense\t$antisense\t$undetermined" ]
                }
                .set { ch_fail_strand }

            MULTIQC_CUSTOM_STRAND_CHECK (
                ch_fail_strand.collect()
            )
            .set { ch_fail_strand_multiqc }
        }
    }

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowRnaseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_fail_mapping_multiqc.ifEmpty([]),
            ch_fail_strand_multiqc.ifEmpty([]),
            FASTQC_TRIM.out.zip.collect{it[1]}.ifEmpty([]), 
            FASTQC_UMITOOLS_TRIMMOMATIC.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            ch_star_multiqc.collect{it[1]}.ifEmpty([]),
            ch_qualimap_multiqc.collect{it[1]}.ifEmpty([]),
            ch_dupradar_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bamstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_inferexperiment_multiqc.collect{it[1]}.ifEmpty([]),
            ch_innerdistance_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionannotation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionsaturation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readdistribution_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readduplication_multiqc.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, fail_percent_mapped)
    }
    NfcoreTemplate.summary(workflow, params, log, fail_percent_mapped, pass_percent_mapped)
}

/*
========================================================================================
    THE END
========================================================================================
*/

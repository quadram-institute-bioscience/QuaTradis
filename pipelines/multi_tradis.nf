#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
      nextflow run multi_tradis.nf --fastq_list fastqs.txt
    Mandatory arguments:
      --fastqs [file]                 A file containing a list of fastqs for processing (each fastq can be gzipped or uncompressed)
      --reference [file]              The reference containing sequences in fasta format
      --refname [str]                 The name for the reference
      --readlen [int]                 The average read length across the experiment
      --aligner [str]                 The alignment / mapping tool to use (bwa, smalt, minimap2)
      --threads [int]                 The number of threads to use per tradis instance (default: 1)
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, awsbatch, test
    Generic
      --outdir [file]                 The output directory where the results will be saved (Default: './results')
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic (Default: false)
      --help                          Show this message
    """.stripIndent()
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Make sure output directory is specified as a file
outdir = file(params.outdir)
STATS_FILE = file("${outdir}/quatradis.stats")


// Check inputs and make channels
if (params.reference) { ch_reference = file(params.reference, checkIfExists: true) } else { exit 1, 'Reference FastA file not specified, or not found on filesystem.' }
if (params.fastqs) { fastq_list = file(params.fastqs, checkIfExists: true).readLines(); ch_fastqs = Channel.fromPath( fastq_list ) } else { exit 1, 'Fastq list file not specified, or not found on filesystem.' }

// Setup command line options if they are required
opt_aligner = ""
opt_threads = ""
opt_tag = ""
opt_mismatch = ""
opt_mapping_score = ""
if (params.aligner) { opt_aligner = "--aligner=${params.aligner}" }
if (params.threads) { opt_threads = "--threads=${params.threads}" }
if (params.tag) { opt_tag = "--tag=${params.tag}" }
if (params.mismatch) { opt_mismatch = "--mismatch=${params.mismatch}" }
if (params.mapping_score) { opt_mapping_score = "--tag=${params.mapping_score}" }



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       HEADER LOG INFO                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

log.info "-\033[2m--------------------------------------------------\033[0m-"
log.info "Configuration Summary"
def summary = [:]
summary['Run Name']               = custom_runName ?: workflow.runName
summary['Output directory']       = params.outdir
summary['Fastq list file']        = params.fastqs
summary['Indexed reference file'] = params.reference
summary['Reference name']         = params.refname
summary['Read Length']            = params.readlen
summary['Tag']                    = params.tag
summary['Aligner']                = params.aligner
summary['Threads per tradis']     = params.threads
summary['Config files']           = workflow.configFiles.join(', ')
log.info summary.collect { k,v -> "- ${k.padRight(20)}: $v" }.join('\n')
log.info "-\033[2m--------------------------------------------------\033[0m-"




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       START THE PIPELINE                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


log.info "Tradis nextflow pipeline starting\n"

// Get sequence length for contig
process indexReference {
    publishDir "${params.outdir}/ref/", mode: "rellink", overwrite: true

    input:
    file reference_fa from ch_reference

    output:
    file "*"
    val 'ok' into ch_index_reference

    """
    index_reference ${reference_fa} ${params.refname}
    """
}

process tradis {
    publishDir "${params.outdir}/tradis/${fq}", mode: "rellink", overwrite: true

    input:
    file fq from ch_fastqs
    file reference_fa from ch_reference
    val 'ok' from ch_index_reference

    output:
    file "tradis_out.plot.stats" into ch_stats1
    file "tradis_out.plot.stats" into ch_stats2

    file "*"

    """
    tradis ${fq} ${reference_fa} --output_prefix=tradis_out ${opt_aligner} ${opt_threads} ${opt_tag} ${opt_mismatch} ${opt_mapping_score}
    """
}

// Collects all line count checks and blends them into a single file
ch_stats1
    .collectFile(name: "stats_header", newLine: true) { statsfile -> statsfile.text.split('\n')[0] }
    .set { ch_merged_stats1 }

ch_stats2
    .collectFile(name: "stats_content", newLine: true) { statsfile -> statsfile.text.split('\n')[1] }
    .set { ch_merged_stats2 }


process combineStats {

    input:
    file stats_file1 from ch_merged_stats1
    file stats_file2 from ch_merged_stats2

    """
    head -n 1 ${stats_file1} > ${STATS_FILE} && cat ${stats_file2} >> ${STATS_FILE}
    """
}
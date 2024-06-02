/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowCrisprseq.initialise(params, log)
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo )   : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

//
// MODULE
//
include { FIND_ADAPTERS         } from '../modules/local/find_adapters'
include { ORIENT_REFERENCE      } from '../modules/local/orient_reference'
include { CIGAR_PARSER          } from '../modules/local/cigar_parser'
include { PREPROCESSING_SUMMARY } from '../modules/local/preprocessing_summary'
include { CLUSTERING_SUMMARY    } from '../modules/local/clustering_summary'
include { ALIGNMENT_SUMMARY     } from '../modules/local/alignment_summary'
include { TEMPLATE_REFERENCE    } from '../modules/local/template_reference'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                   } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS               } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { PEAR                                      } from '../modules/nf-core/pear/main'
include { CAT_FASTQ                                 } from '../modules/nf-core/cat/fastq/main'
include { SEQTK_SEQ as SEQTK_SEQ_MASK               } from '../modules/nf-core/seqtk/seq/main'
include { SEQTK_SEQ as SEQTK_SEQ_FATOFQ             } from '../modules/nf-core/seqtk/seq/main'
include { VSEARCH_CLUSTER                           } from '../modules/nf-core/vsearch/cluster/main'
include { VSEARCH_SORT                              } from '../modules/nf-core/vsearch/sort/main'
include { RACON as RACON_1                          } from '../modules/nf-core/racon/main'
include { RACON as RACON_2                          } from '../modules/nf-core/racon/main'
include { BOWTIE2_ALIGN                             } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD                             } from '../modules/nf-core/bowtie2/build/main'
include { BWA_MEM                                   } from '../modules/nf-core/bwa/mem/main'
include { BWA_INDEX                                 } from '../modules/nf-core/bwa/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ORIGINAL } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_UMI_1    } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_UMI_2    } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TEMPLATE } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX                            } from '../modules/nf-core/minimap2/index/main'
include { MEDAKA                                    } from '../modules/nf-core/medaka/main'
include { CUTADAPT                                  } from '../modules/nf-core/cutadapt/main'
include { SAMTOOLS_INDEX                            } from '../modules/nf-core/samtools/index/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE GROOVY FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CRISPRSEQ_FGC {

    ch_versions = Channel.empty()


    //
    // Create input channel from input file provided through params.input
    //

    //
    // Create input channel from input file provided through params.input
    //
    Channel.fromSamplesheet("input")
    .multiMap { meta, fastq_1, fastq_2, x, y, z ->
        // x (reference), y (protospacer), and z (template) are part of the targeted workflows and we don't need them
        reads:   [ meta, fastq_1 ]
    }
    .set { ch_input }

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_input.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    ch_input_cutadapt = ch_input.combine(Channel.value([[]]))

    ch_trimmed = Channel.empty()

    if (params.overrepresented) {
        //
        // MODULE: Find overrepresented sequences
        FIND_ADAPTERS (
            FASTQC.out.zip
        )
        .adapters
        .join(ch_input)
        .groupTuple(by: [0])
        // Separate samples by containing overrepresented sequences or not
        .branch {
            meta, adapter_seqs, reads ->
                no_adapters: adapter_seqs[0].size() == 0
                    return [ meta, reads[0] ]
                adapters   : adapter_seqs[0].size() > 0
                    return [ meta, reads[0], adapter_seqs[0] ]
        }
        .set { ch_adapter_seqs }
        ch_versions = ch_versions.mix(FIND_ADAPTERS.out.versions.first())


        //
        // MODULE: Trim adapter sequences
        //
        CUTADAPT (
            ch_adapter_seqs.adapters
        )
        ch_versions = ch_versions.mix(CUTADAPT.out.versions)

        ch_adapter_seqs.no_adapters
        .mix(CUTADAPT.out.reads)
        .groupTuple(by: [0])
        .map {
            meta, fastq ->
                    return [ meta, fastq.flatten() ]
        }
        .set{ ch_trimmed }
    } else {
        ch_trimmed = ch_input
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

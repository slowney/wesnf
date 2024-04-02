/*
-----------------------------------------
    PRINT PARAMS SUMMARY
-----------------------------------------
*/
include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowWESnf.initialise(params, log)

/*
-----------------------------------------
    CONFIG FILES
-----------------------------------------
*/
ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
-----------------------------------------
    IMPORT MODULES
-----------------------------------------
*/
include { FASTQC                      }                 from '../modules/nf-core/fastqc/main'
include { MULTIQC                     }                 from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS }                 from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SEQTK_SAMPLE }                                from '../modules/nf-core/seqtk/sample/main'
include { BWAMEM2_MEM }                                 from '../modules/nf-core/bwamem2/mem/main'
include { PICARD_SORTSAM }                              from '../modules/nf-core/picard/sortsam/main'
include { PICARD_MARKDUPLICATES }                       from '../modules/nf-core/picard/markduplicates/main'
include { PICARD_BEDTOINTERVALLIST as TargetInterval}   from '../modules/nf-core/picard/bedtointervallist/main'
include { VARDICTJAVA }                                 from '../modules/nf-core/vardictjava/main'

/*
-----------------------------------------
    IMPORT SUBWORKFLOWS
-----------------------------------------
*/
include { INPUT_CHECK }                 from '../subworkflows/local/input_check'
include { BAM_QC_PICARD }               from '../subworkflows/nf-core/bam_qc_picard/main'
include { PrepGenome }                  from '../subworkflows/local/PrepGenome/main'

/*
-----------------------------------------
    RUN MAIN WORKFLOW
-----------------------------------------
*/

fasta       = params.fasta
                ? Channel.fromPath(params.fasta).first()
                : Channel.empty()

fasta_fai   = params.fasta_fai                  
                ? Channel.fromPath(params.fasta_fai).collect()                
                : Channel.empty()

fasta_dict  = params.fasta_dict                 
                ? Channel.fromPath(params.fasta_dict).collect()           
                : Channel.empty()

target_bed  = params.target_bed
                ? Channel.fromPath(params.target_bed).collect()
                : Channel.empty()

// Info required for completion email and summary
def multiqc_report = []

workflow WESNF {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowWESnf.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowWESnf.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

    // Configure inputs
    fqIn_ch = Channel.fromSamplesheet("input").map{ meta, read1, read2 -> [meta, [read1, read2]] }

    PrepGenome(fasta, fasta_fai)

    // Gather any new genome files from PrepGenome that aren't specified in params
    bwamem2                 = params.bwamem2                ? Channel.fromPath(params.bwamem2).collect()
                                                            : PrepGenome.out.bwamem2
    fasta_fai               = params.fasta_fai              ? Channel.fromPath(params.fasta_fai).collect()
                                                            : PrepGenome.out.fasta_fai
    fasta_dict              = params.fasta_dict             ? Channel.fromPath(params.fasta_dict).collect()
                                                            : PrepGenome.out.fasta_dict

    // Downsample input fqs
    DS_ch = fqIn_ch.branch{
                Sample: it[0]["DownsampleDepth"] && it[0]["DownsampleDepth"] != null
                No_Sample: it[0]["DownsampleDepth"] == null | !("DownsampleDepth" in it[0]) }

    SEQTK_SAMPLE(
        DS_ch.Sample
        .map{ meta, reads -> [meta, reads, meta.DownsampleDepth] } )

    wfIn_ch = DS_ch.No_Sample.mix(SEQTK_SAMPLE.out.reads)

    // Align and Mark Dup - MarkDup needs coordinate sorted bam
    BWAMEM2_MEM(wfIn_ch, bwamem2.map{ it -> [ [ id:'index' ], it ] }, false)

    PICARD_SORTSAM(BWAMEM2_MEM.out.bam, 'coordinate')

    PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam, fasta.map{ it -> [[ id:'fasta' ], it ]}, fasta_fai.map{ it -> [[id:'fai'], it ]})

    // Build interval files for picard metrics
    TargetInterval(
            target_bed.map{target_bed -> [[ id:target_bed.baseName ], target_bed]},
            fasta_dict.map{ it -> [[ id:'dict' ], it ]},
            [])


    BAM_QC_in = PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai)

    // HS Metrics is slightly adjusted to take in the target intervals twice so we don't
    // have to use a probe bed as input.
    BAM_QC_in = BAM_QC_in.map{
        meta, bam, bai ->
        [meta, bam, bai, TargetInterval.out.interval_list.get()[1], TargetInterval.out.interval_list.get()[1]]
    }

    BAM_QC_PICARD(BAM_QC_in, fasta.map{ it -> [[ id:'fasta' ], it ]}, fasta_fai.map{ it -> [[id:'fai'], it ]}, fasta_dict.map{ it -> [[ id:'dict' ], it ]})

    // // Variant Calling
    varCallIn_ch = PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai)
                        .map{meta, bam, bai -> [meta, bam, bai, target_bed.get()]}

    VARDICTJAVA(varCallIn_ch, fasta.map{ it -> [[ id:'fasta' ], it ]}, fasta_fai.map{ it -> [[id:'fai'], it ]})
}

/*
-----------------------------------------
    COMPLETION EMAIL AND SUMMARY
-----------------------------------------
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
-----------------------------------------
    THE END
-----------------------------------------
*/

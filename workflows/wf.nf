/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_metafast_pipeline'

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { MEGAHIT                } from '../modules/nf-core/megahit/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { METABAT2_METABAT2      } from '../modules/nf-core/metabat2/metabat2/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METAFAST {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    /*
    ================================================================================
                                    QC Reports for reads
    ================================================================================
    */

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run Fastp for trimming and filtering
    //
    FASTP (
        ch_samplesheet,
        params.adapter_fasta ?: [],
        params.discard_trimmed_pass ?: false,
        params.save_trimmed_fail ?: false,
        params.save_merged ?: false
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    /*
    ================================================================================
                                    Assembly and Mapping
    ================================================================================
    */

    //
    // MODULE: Run MEGAHIT for assembly
    //
    MEGAHIT (
        FASTP.out.reads
    )
    ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
/*
    //
    // MODULE: Index the assembled contigs with BWA
    //
    BWA_INDEX (
        MEGAHIT.out.contigs
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    //
    // MODULE: Align reads to contigs with BWA MEM
    //
    BWA_MEM (
        FASTP.out.reads,
        BWA_INDEX.out.index
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    //
    // MODULE: Sort BAM files with SAMtools
    //
    SAMTOOLS_SORT (
        BWA_MEM.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    //
    // MODULE: Index sorted BAM files with SAMtools
    //
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    /*
    ================================================================================
                                    Binning
    ================================================================================
    */
/*
    //
    // MODULE: Summarize BAM contig depths
    //
    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS (
        SAMTOOLS_SORT.out.bam,
        SAMTOOLS_INDEX.out.bai
    )
    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first())

    //
    // MODULE: Run MetaBAT2 for binning
    //
    METABAT2_METABAT2 (
        MEGAHIT.out.contigs,
        METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depths
    )
    ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions.first())
*/
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    /*
    ================================================================================
                                    MultiQC-Summary
    ================================================================================
    */

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()

    summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

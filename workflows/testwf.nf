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
                                    Assembly
    ================================================================================
    */

    //
    // MODULE: Run MEGAHIT for assembly
    //
    MEGAHIT (
        FASTP.out.reads
    )
    ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())

    // Prepare assembled contigs for downstream analysis
    ch_assemblies = MEGAHIT.out.contigs.map { meta, contigs ->
        [meta + [assembler: 'megahit'], contigs]
    }

    /*
    ================================================================================
                                    Post-assembly steps
    ================================================================================
    */

    //
    // MODULE: Run BWA indexing on MEGAHIT contigs
    //
    BWA_INDEX (
        ch_assemblies
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    //
    // MODULE: Run BWA alignment of reads to contigs
    //
    BWA_MEM (
        FASTP.out.reads,
        BWA_INDEX.out.index,
        true
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    //
    // MODULE: Sort BAM files
    //
    SAMTOOLS_SORT (
        BWA_MEM.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    //
    // MODULE: Index sorted BAM files
    //
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // MODULE: Summarize BAM contig depths
    //
    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS (
        SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
    )
    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first())

    //
    // MODULE: Run MetaBAT2 for binning
    //
    METABAT2_METABAT2 (
        ch_assemblies,
        METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    )
    ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions.first())

    // ... (continue with any additional steps in your workflow)

    emit:
    assemblies = ch_assemblies
    versions   = ch_versions
    multiqc    = ch_multiqc_files
}

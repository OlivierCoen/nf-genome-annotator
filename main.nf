#!/usr/bin/env nextflow

nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    genomeannotator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/genomeannotator
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOMEANNOTATOR         } from './workflows/genomeannotator'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_genomeannotator_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_genomeannotator_pipeline'

include { getField                } from './subworkflows/local/utils_nfcore_genomeannotator_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_GENOMEANNOTATOR {

    take:
    samplesheet

    main:

    //
    // WORKFLOW: Run pipeline
    //
    GENOMEANNOTATOR( samplesheet )

    emit:
    results = GENOMEANNOTATOR.out.results
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_GENOMEANNOTATOR (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_GENOMEANNOTATOR.out.results.map{ rec -> rec.multiqc_report }.toList()
    )

    publish:
    results            = NFCORE_GENOMEANNOTATOR.out.results
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {

    results {
        path { rec ->
            rec.final_annotation >> "${rec.id}/${rec.id}.annotation.gff3"
            rec.proteome         >> "${rec.id}/${rec.id}.proteome.faa"
            rec.masked_genome    >> "${rec.id}/${rec.id}.masked_genome.fna" 
            rec.multiqc_report   >> "${rec.id}/${rec.id}.multiqc_report.html"

            rec.structural_annotation    >> "${rec.id}/structural_annotations/final/"
            rec.intermediate_annotations >> "${rec.id}/structural_annotations/partially_cleaned/"
            rec.alternative_annotations  >> "${rec.id}/structural_annotations/alternative/"
           

            rec.eggnogmapper_output   >> "${rec.id}/functional_annotations/eggnog_mapper/"
            rec.interproscan_output   >> "${rec.id}/functional_annotations/interproscan/"
            rec.functional_annotation >> "${rec.id}/functional_annotations/"

            rec.genome_stats                >> "${rec.id}/quality_controls/"
            rec.structural_annotation_stats >> "${rec.id}/quality_controls/${rec.id}.structural_annotation_stats.yaml"
            rec.functional_annotation_stats >> "${rec.id}/quality_controls/${rec.id}.functional_annotation_stats.yaml"
            rec.omark                       >> "${rec.id}/quality_controls/"
        }
    }

}

/*
rec.reads.each{ file -> file >> "${rec.id}/structural_annotations/mappings/downloaded_reads/" }
rec.bams.each { file -> file >> "${rec.id}/structural_annotations/mappings/" }
rec.bais.each { file -> file >> "${rec.id}/structural_annotations/mappings/" }
*/
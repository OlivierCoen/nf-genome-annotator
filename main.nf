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
    final_annotation            = GENOMEANNOTATOR.out.final_annotation
    genome_stats                = GENOMEANNOTATOR.out.genome_stats
    masked_genome               = GENOMEANNOTATOR.out.masked_genome
    downloaded_reads            = GENOMEANNOTATOR.out.downloaded_reads
    mapping                     = GENOMEANNOTATOR.out.mapping
    structural_annotation       = GENOMEANNOTATOR.out.structural_annotation
    intermediate_annotations    = GENOMEANNOTATOR.out.intermediate_annotations
    alternative_annotations     = GENOMEANNOTATOR.out.alternative_annotations
    proteome                    = GENOMEANNOTATOR.out.proteome
    eggnogmapper_output         = GENOMEANNOTATOR.out.eggnogmapper_output
    interproscan_output         = GENOMEANNOTATOR.out.interproscan_output
    functional_annotation       = GENOMEANNOTATOR.out.functional_annotation
    structural_annotation_stats = GENOMEANNOTATOR.out.structural_annotation_stats
    functional_annotation_stats = GENOMEANNOTATOR.out.functional_annotation_stats
    omark_results               = GENOMEANNOTATOR.out.omark_results
    multiqc_report              = GENOMEANNOTATOR.out.multiqc_report
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
        NFCORE_GENOMEANNOTATOR.out.multiqc_report
    )

    publish:
    final_annotation            = NFCORE_GENOMEANNOTATOR.out.final_annotation
    genome_stats                = NFCORE_GENOMEANNOTATOR.out.genome_stats
    masked_genome               = NFCORE_GENOMEANNOTATOR.out.masked_genome
    downloaded_reads            = NFCORE_GENOMEANNOTATOR.out.downloaded_reads
    mapping                     = NFCORE_GENOMEANNOTATOR.out.mapping
    structural_annotation       = NFCORE_GENOMEANNOTATOR.out.structural_annotation
    intermediate_annotations    = NFCORE_GENOMEANNOTATOR.out.intermediate_annotations
    alternative_annotations     = NFCORE_GENOMEANNOTATOR.out.alternative_annotations
    proteome                    = NFCORE_GENOMEANNOTATOR.out.proteome
    eggnogmapper_output         = NFCORE_GENOMEANNOTATOR.out.eggnogmapper_output
    interproscan_output         = NFCORE_GENOMEANNOTATOR.out.interproscan_output
    functional_annotation       = NFCORE_GENOMEANNOTATOR.out.functional_annotation
    structural_annotation_stats = NFCORE_GENOMEANNOTATOR.out.structural_annotation_stats
    functional_annotation_stats = NFCORE_GENOMEANNOTATOR.out.functional_annotation_stats
    omark_results               = NFCORE_GENOMEANNOTATOR.out.omark_results
    multiqc_report              = NFCORE_GENOMEANNOTATOR.out.multiqc_report

    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {

    final_annotation {
        path { meta, file ->
            file >> "${meta.id}/final/${meta.id}.annotation.gff3"
        }
    }

    proteome {
        path { meta, file ->
            file >> "${meta.id}/final/${meta.id}.proteome.faa"
        }
    }

    masked_genome {
        path { meta, file ->
            file >> "${meta.id}/final/${meta.id}.masked_genome.fna"
        }
    }


    
    downloaded_reads {
        path { meta, file ->
            file >> "${meta.id}/structural_annotations/mappings/downloaded_reads/"
        }
    }

    mapping {
        path { meta, bam, bai ->
            bam >> "${meta.id}/structural_annotations/mappings/"
            bai >> "${meta.id}/structural_annotations/mappings/"
        }
    }

    structural_annotation {
        path { meta, file ->
            file >> "${meta.id}/structural_annotations/final/"
        }
    }

    intermediate_annotations {
        path { meta, file ->
            file >> "${meta.id}/structural_annotations/partially_cleaned/"
        }
    }

    alternative_annotations {
        path { meta, file ->
            file >> "${meta.id}/structural_annotations/alternative/"
        }
    }

    

    eggnogmapper_output {
        path { meta, file ->
            file >> "${meta.id}/functional_annotations/eggnog_mapper/"
        }
    }

    interproscan_output {
        path { meta, file ->
            file >> "${meta.id}/functional_annotations/interproscan/"
        }
    }

    functional_annotation {
        path { meta, file ->
            file >> "${meta.id}/functional_annotations/"
        }
    }


    

    genome_stats {
        path { meta, file ->
            file >> "${meta.id}/quality_controls/${meta.id}.genome_stats.${file.extension}"
        }
    }


    structural_annotation_stats {
        path { meta, file ->
            file >> "${meta.id}/quality_controls/${meta.id}.structural_annotation_stats.${file.extension}"
        }
    }

    functional_annotation_stats {
        path { meta, file ->
            file >> "${meta.id}/quality_controls/${meta.id}.functional_annotation_stats.${file.extension}"
        }
    }

    omark_results {
        path { meta, file ->
            file >> "${meta.id}/quality_controls/"
        }
    }

    multiqc_report {
        path { file ->
            file >> "reporting/"
        }
    }

}

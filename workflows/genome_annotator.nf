/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AGAT_SPCOMPLEMENTANNOTATIONS as COMPLEMENT_ANNOTATIONS        } from '../modules/local/agat/spcomplementannotations'
include { AGAT_SPEXTRACTSEQUENCES as GET_PROTEOME                       } from '../modules/local/agat/spextractsequences'

include { GENOME_PREPARATION                                            } from '../subworkflows/local/genome_preparation'
include { GENOME_MASKING                                                } from '../subworkflows/local/genome_masking'
include { STRUCTURAL_ANNOTATION                                         } from '../subworkflows/local/structural_annotation'
include { CLEAN_ANNOTATIONS                                             } from '../subworkflows/local/clean_annotations'
include { FUNCTIONAL_ANNOTATION                                         } from '../subworkflows/local/functional_annotation'
include { QUALITY_CONTROLS                                              } from '../subworkflows/local/qc'
include { MULTIQC_WORKFLOW                                              } from '../subworkflows/local/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOME_ANNOTATOR {

    take:
    ch_genome // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GENOME PREPARATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GENOME_PREPARATION ( ch_genome )
    ch_genome = GENOME_PREPARATION.out.prepared_genome

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GENOME MASKING
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !params.skip_masking ) {
        GENOME_MASKING ( ch_genome )
        ch_genome = GENOME_MASKING.out.masked_genome
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STRUCTURAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    STRUCTURAL_ANNOTATION ( ch_genome )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // COMPLEMENTATION OF ANNOTATION ()WHEN NECESSARY)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_branched_annotations = STRUCTURAL_ANNOTATION.out.annotations
                                .branch{
                                    meta, annotation ->
                                        to_complement: meta.ref_gff != []
                                            [ meta, meta.ref_gff, annotation ]
                                        leave_me_alone: meta.ref_gff == []
                                            [ meta, annotation ]
                                }

    COMPLEMENT_ANNOTATIONS ( ch_branched_annotations.to_complement, [] )

    ch_annotation = ch_branched_annotations.leave_me_alone
                        .mix( COMPLEMENT_ANNOTATIONS.out.gff )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CLEANING OF GTF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    CLEAN_ANNOTATIONS (
        ch_annotation,
        ch_genome
     )
     ch_gff = CLEAN_ANNOTATIONS.out.gff
                .map {
                    meta, file ->
                        [ meta + [main_annotation: true], file ]
                }
     ch_cleaned_gffs = CLEAN_ANNOTATIONS.out.cleaned_gffs

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE PROTEOME
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GET_PROTEOME (
        ch_cleaned_gffs.join( ch_genome ),
        params.codon_usage_id,
        []
    )
    ch_proteome = GET_PROTEOME.out.proteins
                    .filter {
                        meta, file ->
                            meta.main_annotation == true
                    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FUNCTIONAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !params.skip_functional_annotation ) {
        FUNCTIONAL_ANNOTATION (
            ch_proteome,
            ch_gff
        )
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // QUALITY CONTROLS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    QUALITY_CONTROLS (
        ch_genome,
        ch_gff,
        ch_proteome
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MULTIQC
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_multiqc_files = ch_multiqc_files
                        .mix ( QUALITY_CONTROLS.out.busco_short_summaries )
                        .mix ( QUALITY_CONTROLS.out.gtf_stats )
                        .map { meta, file -> file }


    MULTIQC_WORKFLOW(
        ch_multiqc_files,
        ch_versions
    )


    emit:
    multiqc_report = MULTIQC_WORKFLOW.out.report.toList()
    versions       = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

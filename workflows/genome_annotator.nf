/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AGAT_SPCOMPLEMENTANNOTATIONS as COMPLEMENT_ANNOTATIONS        } from '../modules/local/agat/spcomplementannotations'
include { AGAT_SPEXTRACTSEQUENCES as GET_PROTEOME                       } from '../modules/local/agat/spextractsequences'

include { GENOME_MASKING                                                } from '../subworkflows/local/genome_masking'
include { STRUCTURAL_ANNOTATION                                         } from '../subworkflows/local/structural_annotation'
include { CLEAN_ANNOTATIONS                                             } from '../subworkflows/local/clean_annotations'
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

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GENOME MASKING
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !params.skip_masking ) {
        GENOME_MASKING ( ch_genome )
        GENOME_MASKING.out.masked_genome.set { ch_genome }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STRUCTURAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    STRUCTURAL_ANNOTATION ( ch_genome )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // COMPLEMENTATION OF ANNOTATION ()WHEN NECESSARY)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    STRUCTURAL_ANNOTATION.out.gtf
        .branch{
            meta, gtf ->
                to_complement: meta.ref_gff != null
                    [ meta, meta.ref_gff, gtf ]
                leave_me_alone: ref_gff == null
                    [ meta, gtf ]
        }
        .set { ch_branched_gtf }

    COMPLEMENT_ANNOTATIONS ( ch_branched_gtf.to_complement, [] )

    ch_branched_gtf.leave_me_alone
        .mix( COMPLEMENT_ANNOTATIONS.out.gff )
        .set { ch_all_gtfs }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CLEANING OF GTF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    CLEAN_ANNOTATIONS (
        ch_all_gtfs,
        ch_genome
     )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE PROTEOME
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GET_PROTEOME (
        CLEAN_ANNOTATIONS.out.gtf.join ( ch_genome ),
         params.codon_usage_id,
         []
    )


    ch_versions
        .mix ( GENOME_MASKING.out.versions )
        .mix ( STRUCTURAL_ANNOTATION.out.versions )
        .mix ( CLEAN_ANNOTATIONS.out.versions )
        .set { ch_versions }

    MULTIQC_WORKFLOW( ch_versions )


    emit:
    multiqc_report = MULTIQC_WORKFLOW.out.report.toList()
    versions       = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BUSCO_DOWNLOADPROTEINS as DOWNLOAD_ORTHODB_PROTEINS   } from '../../../modules/local/busco/download'
include { BRAKER3                                               } from '../../../modules/local/braker3'
include { TSEBRA_TSEBRA          as TSEBRA                      } from '../../../modules/local/tsebra/tsebra'

include { MMSEQS_DATABASES                                      } from '../../../modules/nf-core/mmseqs/databases'
include { SAMTOOLS_MERGE                                        } from '../../../modules/nf-core/samtools/merge'
include { METAEUK_EASYPREDICT                                   } from '../../../modules/local/metaeuk/easypredict'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow STRUCTURAL_ANNOTATION {

    take:
    ch_genome
    ch_proteins
    ch_bam
    ch_gtf
    ch_hintsfile
    structural_annotator
    species

    main:

    ch_versions = channel.empty()

    if ( structural_annotator == "braker3" ) {

        // ----------------------------------------------------------
        // MERGE MULTIPLE BAM FILES
        // ----------------------------------------------------------

        ch_branched_bam = ch_bam
                            .branch{
                                meta, bams ->
                                    leave_me_alone: bams.size() == 1
                                        [ meta, bams[0] ]
                                    merge_me: bams.size() > 1
                                        [ meta, bams ]
                            }

        SAMTOOLS_MERGE(
            ch_branched_bam.merge_me,
            []
        )

        ch_single_bam = ch_branched_bam.leave_me_alone
                            .mix( SAMTOOLS_MERGE.out.bam )

        // ----------------------------------------------------------
        // RUN BRAKER3 & TSEBRA
        // ----------------------------------------------------------

        ch_braker_input = ch_genome
                            .join( ch_proteins )
                            .join( ch_single_bam )

        def species_arg = species ?: []
        BRAKER3(
            ch_braker_input,
            species_arg
        )

        // MERGE ANNOTATIONS WHEN NECESSARY

        // separate inputs that need to be merged from the rest
        ch_braker_gtfs = BRAKER3.out.gtf
                            .join ( BRAKER3.out.hintsfile )
                            .join( ch_gtf )
                            .join( ch_hintsfile )
                            .map {
                                meta, braker_gtf, braker_hintsfile, gtf, hintsfile ->
                                    if ( gtf != [] && hintsfile == [] ) {
                                        log.warn("Cannot merge BRAKER output with existing gtf for sample ${meta.id} because hintsfile is absent")
                                    } else if ( gtf == [] && hintsfile != [] ) {
                                        log.warn("Cannot merge BRAKER output with existing gtf for sample ${meta.id} because gtf is absent")
                                    }
                            }
                            .branch {
                                meta, braker_gtf, braker_hintsfile, gtf, hintsfile ->
                                    to_merge: gtf != [] && hintsfile != []
                                        [ meta, [ braker_gtf, gtf ], [ braker_hintsfile, hintsfile ] ]
                                    not_to_merge: gtf == [] || hintsfile == []
                                        [ meta, braker_gtf ]
                            }

        // TSEBRA
        TSEBRA(
            ch_braker_gtfs.to_merge,
            [],
            []
        )

        ch_annotations = ch_braker_gtfs.not_to_merge
                                    .mix( TSEBRA.out.merged_gtf )

    } else {

        // ----------------------------------------------------------
        // RUN METAEUK
        // ----------------------------------------------------------

        METAEUK_EASYPREDICT( ch_to_annotate )
        ch_annotations = METAEUK_EASYPREDICT.out.gff

    }

    emit:
    annotations             = ch_annotations
    versions                = ch_versions
}

include { BRAKER3                                               } from '../../../modules/local/braker3'
include { TSEBRA_TSEBRA as TSEBRA                               } from '../../../modules/local/tsebra/tsebra'
include { SAMTOOLS_INDEX                                        } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_MERGE                                        } from '../../../modules/local/samtools/merge'
include { METAEUK_EASYPREDICT                                   } from '../../../modules/local/metaeuk/easypredict'

include { PROTEIN_DB_PREPARATION                                } from '../protein_db_preparation'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow STRUCTURAL_ANNOTATION {

    take:
    ch_genome
    ch_proteins
    ch_grouped_bam_bai
    ch_gtf
    ch_hintsfile
    structural_annotator
    species
    clade
    excluded_clades
    excluded_species

    main:

    ch_versions = channel.empty()

    if ( structural_annotator == "braker3" ) {

        // ----------------------------------------------------------
        // MERGE MULTIPLE BAM FILES
        // ----------------------------------------------------------

        ch_branched_bam = ch_grouped_bam_bai
                            .branch{
                                meta, bams, bais ->
                                    leave_me_alone: bams.size() == 1
                                        [ meta, bams[0], bais[0] ]
                                    merge_me: bams.size() > 1
                                        [ meta, bams, bais ]
                            }

        SAMTOOLS_MERGE( ch_branched_bam.merge_me )

        ch_single_bam = ch_branched_bam.leave_me_alone
                            .mix( SAMTOOLS_MERGE.out.bam )
                            .map{ meta, bam, bai -> [ meta, bam ] }

        // ----------------------------------------------------------
        // PREPARE PROTEIN DB FROM CLADE-SPECIFIC ORTHODB AND CUSTOM PROTEIN FASTA FILES
        // ----------------------------------------------------------

        PROTEIN_DB_PREPARATION(
            ch_proteins,
            clade,
            excluded_clades,
            excluded_species
        )
        ch_prepared_proteins = PROTEIN_DB_PREPARATION.out.proteins

        // ----------------------------------------------------------
        // RUN BRAKER3
        // ----------------------------------------------------------

        ch_braker_input = ch_genome
                            .map{ meta, genome -> [ [id: meta.id], genome ] }
                            .join( ch_prepared_proteins, remainder: true )
                            .join( ch_single_bam, remainder: true )
                            .map{
                                meta, genome, prot, bam ->
                                    [ meta, genome, prot?: [], bam?: [] ]
                            }

        def species_arg = species ?: []
        BRAKER3(
            ch_braker_input,
            species_arg
        )

        // ----------------------------------------------------------
        // MERGE ANNOTATIONS WHEN NECESSARY
        // ----------------------------------------------------------

        // separate inputs that need to be merged from the rest
        ch_branched_braker_gtfs = BRAKER3.out.gtf
                                    .join ( BRAKER3.out.hintsfile )
                                    .join( ch_gtf, remainder: true )
                                    .join( ch_hintsfile, remainder: true )
                                    .subscribe {
                                        meta, braker_gtf, braker_hintsfile, gtf, hintsfile ->
                                            if ( gtf != null && hintsfile == null ) {
                                                log.warn("Cannot merge BRAKER output with existing gtf for sample ${meta.id} because hintsfile is absent")
                                            } else if ( gtf == null && hintsfile != null ) {
                                                log.warn("Cannot merge BRAKER output with existing gtf for sample ${meta.id} because gtf is absent")
                                            }
                                    }
                                    .branch {
                                        meta, braker_gtf, braker_hintsfile, gtf, hintsfile ->
                                            to_merge: gtf != null && hintsfile != null
                                                [ meta, [ braker_gtf, gtf ], [ braker_hintsfile, hintsfile ] ]
                                            not_to_merge: gtf == null || hintsfile == null
                                                [ meta, braker_gtf ]
                                    }

        // TSEBRA
        TSEBRA(
            ch_branched_braker_gtfs.to_merge,
            [],
            []
        )

        ch_annotations = ch_branched_braker_gtfs.not_to_merge
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

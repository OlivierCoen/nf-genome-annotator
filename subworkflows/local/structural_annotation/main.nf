include { BUSCO_DOWNLOADPROTEINS as DOWNLOAD_ORTHODB_PROTEINS   } from '../../../modules/local/busco/download'
include { BRAKER3                                               } from '../../../modules/local/braker3'
include { TSEBRA_TSEBRA          as TSEBRA                      } from '../../../modules/local/tsebra/tsebra'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow STRUCTURAL_ANNOTATION {

    take:
    ch_genome

    main:

    // ----------------------------------------------------------
    // PREPARE INPUTS
    // ----------------------------------------------------------

    ch_genome
        .branch {
            meta, genome ->
                no_proteins: meta.proteins == []
                    [ meta, genome, meta.bam ]
                proteins: meta.proteins != []
                    [ meta, genome, meta.bam, meta.proteins ]
        }
        .set { ch_prepared }

    // ----------------------------------------------------------
    // DOWNLOAD ORTHODB PROTEINS IF NECESSARY
    // ----------------------------------------------------------

    // trick to exectute the process DOWNLOAD_ORTHODB_PROTEINS
    // only if there are elements in the channel ch_prepared.no_proteins
    def orthodb_lineage = params.braker_orthodb_lineage ?: params.busco_lineage
    ch_prepared.no_proteins
        .take (1)
        .map { orthodb_lineage }
        | DOWNLOAD_ORTHODB_PROTEINS

    // add download proteins to the items without proteins
    ch_prepared.no_proteins
        .combine( DOWNLOAD_ORTHODB_PROTEINS.out.proteins )
        .map {
            meta, genome, bam, download_proteins ->
                [ meta, genome, bam, download_proteins ]
        }
        .set { ch_prepared_with_downloaded_proteins }


    // ----------------------------------------------------------
    // RUN BRAKER3
    // ----------------------------------------------------------

    ch_prepared.proteins
        .mix ( ch_prepared_with_downloaded_proteins )
        .set { ch_braker_input }

    BRAKER3( ch_braker_input )

    // ----------------------------------------------------------
    // MERGE ANNOTATIONS WHEN NECESSARY
    // ----------------------------------------------------------

    BRAKER3.out.gtf
        .branch {
            meta, gtf ->
                to_merge: meta.gtf != [] && meta.hintsfile != []
                not_to_merge: meta.gtf == [] || meta.hintsfile == []
         }
         .set { ch_gtfs }

    ch_gtfs.to_merge
        .join ( BRAKER3.out.hintsfile )
        .map {
            meta, gtf, hintsfile ->
                [ meta, [ gtf, meta.gtf ], [ hintsfile, meta.hintsfile ] ]
        }
        .set { ch_tsebra_input }

    TSEBRA( ch_tsebra_input, [], [] )


    ch_gtfs.not_to_merge
        .mix( TSEBRA.out.merged_gtf )
        .set { ch_all_gtfs }


    emit:
    gtf                     = ch_all_gtfs

}

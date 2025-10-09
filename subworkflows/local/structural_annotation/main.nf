include { BRAKER3                           } from '../../../modules/local/braker3'
include { TSEBRA_TSEBRA as TSEBRA           } from '../../../modules/local/tsebra/tsebra'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow STRUCTURAL_ANNOTATION {

    take:
    ch_genome

    main:

    ch_versions = Channel.empty()

    // ----------------------------------------------------------
    // RUN BRAKER3
    // ----------------------------------------------------------

    BRAKER3( ch_genome, [], [], [], [], [] )

    BRAKER3.out.gtf
        .branch {
            meta, gtf ->
                to_merge: meta.gtf != [] && meta.hintsfile != []
                not_to_merge: meta.gtf == [] || meta.hintsfile == []
         }
         .set { ch_gtfs }

    // ----------------------------------------------------------
    // MERGE ANNOTATIONS WHEN NECESSARY
    // ----------------------------------------------------------

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


    ch_versions
        .mix ( BRAKER3.out.versions )
        .set { ch_versions }


    emit:
    gtf                     = ch_all_gtfs
    versions                = ch_versions

}


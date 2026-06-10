include { MMSEQS_DATABASES                                   } from '../../../modules/nf-core/mmseqs/databases'
include { MMSEQS_CREATEDB                                    } from '../../../modules/local/mmseqs/createdb'
include { MMSEQS_CONCATDBS                                   } from '../../../modules/local/mmseqs/concatdbs'
//include { MMSEQS_EXCLUDE_TOO_SMALL_PROTEINS                  } from '../../../modules/local/mmseqs/exclude_too_small_proteins'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow MMSEQS_DB_PREPARATION {

    take:
    ch_proteins
    mmseqs_db
    skip_mmseqs_db_download
    min_prot_db_seq_length

    main:

    // ----------------------------------------------------------
    // DOWNLOAD MMSEQS DB
    // ----------------------------------------------------------

    if ( !skip_mmseqs_db_download ) {

        MMSEQS_DATABASES( mmseqs_db )
        ch_mmseqs_public_db = MMSEQS_DATABASES.out.database
        
    } else {
        ch_mmseqs_public_db = channel.of( "no_db" )
    }

    // ----------------------------------------------------------
    // CREATE DBS FROM CUSTOM PROTEIN SEQUENCES
    // ----------------------------------------------------------

    MMSEQS_CREATEDB( ch_proteins )
    ch_custom_mmseqs_db = MMSEQS_CREATEDB.out.db

    // ----------------------------------------------------------
    // CONCAT ALL PROTEIN DBS INTO A SINLE ONE
    // ----------------------------------------------------------

    ch_all_db = ch_custom_mmseqs_db
                        .combine( ch_mmseqs_public_db ) // cartesian product: add the public protein db to each item separately
                        .map {
                            meta, custom_db, public_db ->
                                if ( public_data == "no_db" ) {
                                    [ meta, custom_db ]
                                } else {
                                    [ meta, input_fasta_list + [public_data] ]
                                }
                        }

    ch_branched_all_db = ch_all_db
                                .branch {
                                    meta, db_list ->
                                        to_concat: db_list.size() > 1
                                        leave_me_alone: db_list.size() <= 1
                                }

    MMSEQS_CONCATDBS ( ch_branched_all_db.to_concat )

    ch_all_combined_db = ch_branched_all_db.leave_me_alone
                                .mix( MMSEQS_CONCATDBS.out.db )

    // ----------------------------------------------------------
    // CHECK HEADERS OR WHOLE PROTEIN DB AND FIX THEM WHEN NECESSARY
    // ----------------------------------------------------------
    /*
    MMSEQS_EXCLUDE_TOO_SMALL_PROTEINS(
        ch_all_combined_db,
        min_prot_db_seq_length
    )
    */
    

    emit:
    proteins  = MMSEQS_CONCATDBS.out.db
}

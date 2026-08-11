nextflow.enable.types = true

include { MMSEQS_DATABASES                                   } from '../../../modules/local/mmseqs/databases'
include { MMSEQS_CREATEDB                                    } from '../../../modules/local/mmseqs/createdb'
include { MMSEQS_CONCATDBS                                   } from '../../../modules/local/mmseqs/concatdbs'
//include { MMSEQS_EXCLUDE_TOO_SMALL_PROTEINS                  } from '../../../modules/local/mmseqs/exclude_too_small_proteins'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def warnNoProteinsFound(rec: Record){
    log.warn(
        "No proteins found for ${rec.id}. Structural annotation will be skipped for this sample. " +
        "Please provide training proteins for this sample in the samplesheet, or unset the --skip_mmseqs_db_download flag."
    )
}

record Input {
    id: String
    training_proteins: Iterable<Path>
    mmseqs_db: String
}


workflow MMSEQS_DB_PREPARATION {

    take:
    ch_input: Channel<Input>
    mmseqs_db
    skip_mmseqs_db_download
    min_prot_db_seq_length

    main:

    // ----------------------------------------------------------
    // DOWNLOAD MMSEQS DB
    // ----------------------------------------------------------

    if ( !skip_mmseqs_db_download ) {

        // the mmseqs db supplied in the samplesheet takes precedence over the one provided globally in the parameters
        ch_input = ch_input.map{ rec ->
            rec + record(mmseqs_db_to_download: rec.mmseqs_db ?: mmseqs_db)
        }

        MMSEQS_DATABASES(
            ch_input.map { rec -> rec.mmseqs_db_to_download }.unique()
        )

        ch_downloaded_mmseqs_db = MMSEQS_DATABASES.out.map{ rec ->
            record(mmseqs_db_to_download: rec.id, public_mmseqs_db: rec.db)
        }

        ch_input = ch_input.join( ch_downloaded_mmseqs_db, by: 'mmseqs_db_to_download' )

    } else {
        ch_input = ch_input.map { rec -> rec + record(mmseqs_db: null) }
    }

    // ----------------------------------------------------------
    // CREATE DBS FROM CUSTOM PROTEIN SEQUENCES
    // ----------------------------------------------------------

    ch_mmseqs_db_to_create = ch_input.filter { rec -> rec.training_proteins.size() > 0 }
    ch_leave_me_alone      = ch_input.filter { rec -> rec.training_proteins.size() == 0 }

    MMSEQS_CREATEDB(
        ch_mmseqs_db_to_create.map { rec -> record(id: rec.id, sequences: rec.training_proteins) }
    )

    ch_mmseqs_db_created = ch_mmseqs_db_to_create.join( MMSEQS_CREATEDB.out, by: 'id')
    ch_leave_me_alone    = ch_leave_me_alone.map { rec -> rec + record(custom_mmseqs_db: null) }

    ch_input = ch_mmseqs_db_created.mix( ch_leave_me_alone )

    // ----------------------------------------------------------
    // CONCAT ALL PROTEIN DBS INTO A SINLE ONE
    // ----------------------------------------------------------

    ch_input = ch_input.map{ rec ->
        rec + record(mmseqs_dbs: [rec.public_mmseqs_db, rec.custom_mmseqs_db].findAll{ db -> db != null } )
    }

    ch_mmseqs_dbs_to_concat = ch_input.filter { rec -> rec.mmseqs_dbs.size() > 1 }
    ch_single_mmseqs_db     = ch_input.filter { rec -> rec.mmseqs_dbs.size() == 1 }
    ch_no_mmseqs_db         = ch_input.filter { rec -> rec.mmseqs_dbs.size() == 0 }

    MMSEQS_CONCATDBS( ch_mmseqs_dbs_to_concat )

    ch_concatenated_dbs = ch_mmseqs_dbs_to_concat.join( MMSEQS_CONCATDBS.out, by: 'id' )
    ch_single_mmseqs_db = ch_single_mmseqs_db.map { rec -> rec + record(mmseqs_db: rec.mmseqs_dbs[0]) }

    ch_no_mmseqs_db.map{ rec -> warnNoProteinsFound(rec) }

    ch_mmseqs_db = ch_concatenated_dbs.mix( ch_single_mmseqs_db )

    // ----------------------------------------------------------
    // FILTER DB
    // ----------------------------------------------------------
    /*
    MMSEQS_EXCLUDE_TOO_SMALL_PROTEINS(
        ch_all_combined_db,
        min_prot_db_seq_length
    )
    */


    emit:
    db  = ch_mmseqs_db
}

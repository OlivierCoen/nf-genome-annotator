nextflow.enable.types = true

include { MMSEQS_DB_PREPARATION                                 } from '../mmseqs_db_preparation'
include { METAEUK_EASYPREDICT                                   } from '../../../modules/local/metaeuk/easypredict'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    fasta: Path
    training_proteins: Iterable<Path>
    mmseqs_db: String
}


workflow METAEUK {

    take:
    ch_input: Channel<Input>
    mmseqs_db
    skip_mmseqs_db_download
    min_prot_db_seq_length

    main:

    // ----------------------------------------------------------
    // PREPARE MMSEQS PROTEIN DB FROM THE CHOSEN MMSEQS DB AND CUSTOM PROTEIN FASTA FILES
    // ----------------------------------------------------------

    MMSEQS_DB_PREPARATION(
        ch_input.map { rec -> rec.subMap(['id', 'training_proteins', 'mmseqs_db']) },
        mmseqs_db,
        skip_mmseqs_db_download,
        min_prot_db_seq_length
    )

    ch_input = ch_input.join( MMSEQS_DB_PREPARATION.out.db, by: 'id' )

    // ----------------------------------------------------------
    // RUN METAEUK
    // ----------------------------------------------------------

    METAEUK_EASYPREDICT(
        ch_input.map{ rec -> record(id: rec.id, fasta: rec.fasta, db: rec.mmseqs_db)}
    )

    emit:
    annotated = METAEUK_EASYPREDICT.out

}

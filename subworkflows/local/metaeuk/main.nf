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
    training_proteins: Set<Path>
    mmseqs_db: String
}


workflow METAEUK {

    take:
    ch_input: Channel<Input>
    mmseqs_db: String
    skip_mmseqs_db_download: Boolean
    min_prot_db_seq_length: Integer

    main:

    // ----------------------------------------------------------
    // PREPARE MMSEQS PROTEIN DB FROM THE CHOSEN MMSEQS DB AND CUSTOM PROTEIN FASTA FILES
    // ----------------------------------------------------------

    ch_mmseqs_db = MMSEQS_DB_PREPARATION(
        ch_input.map { rec -> rec.subMap(['id', 'training_proteins', 'mmseqs_db']) },
        mmseqs_db,
        skip_mmseqs_db_download,
        min_prot_db_seq_length
    )

    ch_input = ch_input.join( ch_mmseqs_db, by: 'id' )

    // ----------------------------------------------------------
    // RUN METAEUK
    // ----------------------------------------------------------

    ch_metaeuk_out = METAEUK_EASYPREDICT(
        ch_input.map{ rec -> record(id: rec.id, fasta: rec.fasta, db: rec.mmseqs_db)}
    )

    ch_input = ch_input.join(
        ch_metaeuk_out.map { rec -> record(id: rec.id, structural_annotation: rec.gff) },
        by: 'id'
    )

    emit:
    annotated = ch_input

}

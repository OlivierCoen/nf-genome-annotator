include { ORTHODB_MAKECLADEDB                                   } from '../../../modules/local/orthodb/make_clade_db'
include { SEQKIT_CONCAT                                         } from '../../../modules/local/seqkit/concat'
include { CHECK_FASTA                                           } from '../../../modules/local/check_fasta'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow PROTEIN_DB_PREPARATION {

    take:
    ch_proteins
    clade
    excluded_clades
    excluded_species
    skip_orthodb_preparation

    main:

    // ----------------------------------------------------------
    // DOWNLOAD CLASE-SPECIFIC ORTHODB PROTEIN DB
    // ----------------------------------------------------------

    if ( !skip_orthodb_preparation ) {
        ORTHODB_MAKECLADEDB(
            clade,
            excluded_clades,
            excluded_species
        )
        ch_orthodb_proteins = ORTHODB_MAKECLADEDB.out.db
    } else {
        ch_orthodb_proteins = channel.of( [] )
    }

    // ----------------------------------------------------------
    // CONCAT ALL PROTEIN DBS INTO A SINLE ONE
    // ----------------------------------------------------------

    ch_to_concat = ch_proteins
                        .combine( ch_orthodb_proteins ) // cartesian product: add the clade orthodb protein db to each item separately
                        .map {
                            if ( orthodb_data == [] ) {
                                [ meta, input_fasta_list ]
                            } else {
                                [ meta, input_fasta_list + [orthodb_data] ]
                        }

    SEQKIT_CONCAT ( ch_to_concat )
    ch_combined_proteins = SEQKIT_CONCAT.out.fastx

    // ----------------------------------------------------------
    // CHECK HEADERS OR WHOLE PROTEIN DB AND FIX THEM WHEN NECESSARY
    // ----------------------------------------------------------

    CHECK_FASTA(
        SEQKIT_CONCAT.out.fastx,
        type="protein",
        fix_headers=true,
        fix_sequences=true
    )

    emit:
    proteins  = CHECK_FASTA.out.fasta
}

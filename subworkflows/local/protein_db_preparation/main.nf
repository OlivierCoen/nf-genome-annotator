include { ORTHODB_MAKECLADEDB                                   } from '../../../modules/local/orthodb/make_clade_db'
include { SEQKIT_CONCAT                                         } from '../../../modules/local/seqkit/concat'
include { CHECK_FASTA_HEADERS                                   } from '../../../modules/local/check_fasta_headers'


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

    main:

    // ----------------------------------------------------------
    // DOWNLOAD CLASE-SPECIFIC ORTHODB PROTEIN DB
    // ----------------------------------------------------------

    ORTHODB_MAKECLADEDB(
        clade,
        excluded_clades,
        excluded_species
    )

    // ----------------------------------------------------------
    // CONCAT ALL PROTEIN DBS INTO A SINLE ONE
    // ----------------------------------------------------------

    ch_to_concat = ch_proteins
                        .combine( ORTHODB_MAKECLADEDB.out.db ) // cartesian product: add the clade orthodb protein db to each item separately
                        .map {
                            meta, input_fasta_list, orthodb_data ->
                                [ meta, input_fasta_list + [orthodb_data] ]
                        }

    SEQKIT_CONCAT ( ch_to_concat )
    ch_combined_proteins = SEQKIT_CONCAT.out.fastx

    // ----------------------------------------------------------
    // CHECK HEADERS OR WHOLE PROTEIN DB AND FIX THEM WHEN NECESSARY
    // ----------------------------------------------------------

    CHECK_FASTA_HEADERS(
        SEQKIT_CONCAT.out.fastx,
        fix=true
    )

    emit:
    proteins  = CHECK_FASTA_HEADERS.out.fasta
}

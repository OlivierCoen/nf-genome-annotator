include { ORTHODB_MAKECLADEDB                                   } from '../../../modules/local/orthodb/make_clade_db'
include { SEQKIT_CONCAT                                         } from '../../../modules/local/seqkit/concat'
include { CHECK_FASTA                                           } from '../../../modules/local/check_fasta'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow ORTHODB_PREPARATION {

    take:
    ch_proteins
    clade
    excluded_clades
    excluded_species
    skip_orthodb_download
    min_prot_db_seq_length

    main:

    // ----------------------------------------------------------
    // DOWNLOAD CLASE-SPECIFIC ORTHODB PROTEIN DB
    // ----------------------------------------------------------

    if ( !skip_orthodb_download ) {

        ORTHODB_MAKECLADEDB(
            clade,
            excluded_clades,
            excluded_species
        )
        ch_orthodb_proteins = ORTHODB_MAKECLADEDB.out.db
        
    } else {
        ch_orthodb_proteins = channel.of( "no_db" )
    }

    // ----------------------------------------------------------
    // CONCAT ALL PROTEIN DBS INTO A SINLE ONE
    // ----------------------------------------------------------

    ch_all_proteins = ch_proteins
                        .combine( ch_orthodb_proteins ) // cartesian product: add the public protein db to each item separately
                        .map {
                            meta, input_fasta_list, public_data ->
                                if ( public_data == "no_db" ) {
                                    [ meta, input_fasta_list ]
                                } else {
                                    [ meta, input_fasta_list + [public_data] ]
                                }
                        }

    ch_branched_all_proteins = ch_all_proteins
                                .branch {
                                    meta, protein_list ->
                                        to_concat: protein_list.size() > 1
                                        leave_me_alone: protein_list.size() <= 1
                                }

    SEQKIT_CONCAT ( ch_branched_all_proteins.to_concat )

    ch_all_combined_proteins = ch_branched_all_proteins.leave_me_alone
                                .mix( SEQKIT_CONCAT.out.fastx )

    // ----------------------------------------------------------
    // CHECK HEADERS OR WHOLE PROTEIN DB AND FIX THEM WHEN NECESSARY
    // ----------------------------------------------------------

    CHECK_FASTA(
        ch_all_combined_proteins,
        type="protein",
        fix_headers=true,
        fix_sequences=true,
        min_sequence_length=min_prot_db_seq_length
    )

    emit:
    proteins  = CHECK_FASTA.out.fasta
}

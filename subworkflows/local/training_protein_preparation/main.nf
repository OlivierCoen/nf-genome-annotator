nextflow.enable.types = true

include { ORTHODB_MAKECLADEDB                                   } from '../../../modules/local/orthodb/make_clade_db'
include { SEQKIT_CONCAT                                         } from '../../../modules/local/seqkit/concat'
include { CHECK_PROTEIN_FASTA                                   } from '../../../modules/local/check/protein_fasta'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def getOrthodbHash(rec) {
    "${rec.clade}_${rec.excluded_clades.join("_")}_${rec.excluded_species.join("_")}".hashCode()
}

def warnNoProteinsFound(rec: Record){
    log.warn("No proteins found for ${rec.id}. Structural annotation will be skipped for this sample. Please provide training proteins for this sample in the samplesheet, or unset the --skip_orthodb_download flag.")
}

record Input {
    id: String
    clade: String
    orthodb_clade: String
    excluded_clades: Iterable<String>
    excluded_species: Iterable<String>
    training_proteins: Iterable<Path>
}


workflow TRAINING_PROTEIN_PREPARATION {

    take:
    ch_input: Channel<Input>
    skip_orthodb_download
    min_prot_db_seq_length

    main:

    // ----------------------------------------------------------
    // PREPARE ORTHODB PROTEIN DB FROM CLADE-SPECIFIC ORTHODB AND CUSTOM PROTEIN FASTA FILES
    // ----------------------------------------------------------

    if ( !skip_orthodb_download ) {

        ch_input = ch_input
                        .map { rec -> record(
                            clade: rec.orthodb_clade ?: rec.clade,
                            excluded_clades: rec.excluded_clades ?: [],
                            excluded_species: rec.excluded_species ?: []
                        )}

        ORTHODB_MAKECLADEDB(
            ch_input.map { rec -> rec.subMap(['clade', 'excluded_clades', 'excluded_species']) }.unique()
        )

        ch_input = ch_input
                        .map{ rec -> rec + record(orthodb_hash: getOrthodbHash(rec)) }
                        .join(
                            ORTHODB_MAKECLADEDB.out.map { rec -> rec + record(orthodb_hash: getOrthodbHash(rec)) },
                            by: 'orthodb_hash'
                        )
                        .map { rec -> rec.subMap(rec.keySet() - ['orthodb_hash']) }

    } else {
        ch_input = ch_input.map { rec -> rec + record(orthodb_proteins: null) }
    }

    ch_input = ch_input.map { rec ->
        def orthodb_proteins = rec.orthodb_proteins ? [rec.orthodb_proteins] : []
        rec + record(all_training_proteins: rec.training_proteins + orthodb_proteins)
    }

    // ----------------------------------------------------------
    // CONCAT ALL PROTEIN DBS INTO A SINLE ONE
    // ----------------------------------------------------------

    ch_proteins_to_concat = ch_input.filter{ rec -> rec.all_training_proteins.size() > 1 }
    ch_single_protein_db  = ch_input.filter{ rec -> rec.all_training_proteins.size() == 1 }
    ch_no_protein_db      = ch_input.filter{ rec -> rec.all_training_proteins.size() == 0 }

    SEQKIT_CONCAT (
        ch_proteins_to_concat.map { rec -> record(id: rec.id, fasta_files: rec.all_training_proteins) }
    )

    ch_concatenated_proteins = SEQKIT_CONCAT.out
    ch_single_protein_db     = ch_single_protein_db.map { rec -> rec + record(fasta: rec.all_training_proteins[0]) }

    //ch_no_protein_db.map { rec -> warnNoProteinsFound(rec) }
ch_concatenated_proteins.view{ v -> "concat $v"}
ch_single_protein_db.view{ v -> "single $v"}

    //ch_input = ch_concatenated_proteins.mix( ch_single_protein_db )
ch_input = ch_single_protein_db
    // ----------------------------------------------------------
    // CHECK HEADERS OR WHOLE PROTEIN DB AND FIX THEM WHEN NECESSARY
    // ----------------------------------------------------------

    CHECK_PROTEIN_FASTA(
        ch_input,
        min_prot_db_seq_length
    )

    ch_proteins = ch_input
                    .join(
                        CHECK_PROTEIN_FASTA.out.map { rec -> record(id: rec.id, proteins_fasta: rec.fasta) },
                        by: 'id'
                    )


    emit:
    proteins = ch_proteins
}

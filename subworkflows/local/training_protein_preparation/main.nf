nextflow.enable.types = true

include { SEQKIT_CONCAT                                         } from '../../../modules/local/seqkit/concat'
include { CHECK_PROTEIN_FASTA                                   } from '../../../modules/local/check/protein_fasta'

include { ORTHODB_PREPARATION                                   } from '../orthodb_preparation'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


def warnNoProteinsFound(rec: Record){
    log.warn(
        "No proteins found for ${rec.id}. In case you want to add training proteins for better gene evidence, " +
        "please provide training proteins for this sample in the samplesheet, or unset the --skip_orthodb_download flag."
    )
}

record Input {
    id: String
    orthodb_excluded_clades: List<String>
    orthodb_excluded_species: List<String>
    training_proteins: Set<Path>
}


workflow TRAINING_PROTEIN_PREPARATION {

    take:
    ch_input: Channel<Input>
    skip_orthodb_download: Boolean
    orthodb_local: Path?
    min_prot_db_seq_length: Integer

    main:

    // ----------------------------------------------------------
    // PREPARE ORTHODB PROTEIN DB FROM CLADE-SPECIFIC ORTHODB AND CUSTOM PROTEIN FASTA FILES
    // ----------------------------------------------------------

    if ( orthodb_local ) {
    
        ch_orthodb_local = channel.fromPath( orthodb_local, checkExists: true )
                            .map { file -> record(orthodb_proteins: file) }
                            
        ch_input = ch_input.combine( ch_orthodb_local )
                        
    } else {
            
        if ( !skip_orthodb_download ) {
    
            ch_prepared_orthodb = ORTHODB_PREPARATION( 
                ch_input.map { rec -> rec.subMap(['species', 'orthodb_excluded_clades', 'orthodb_excluded_species']) }.unique()
            )
            ch_input = ch_input.join(ch_prepared_orthodb, by: 'species')
    
    
        } else {
            ch_input = ch_input.map { rec -> rec + record(orthodb_proteins: null) }
        }
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

    ch_concatenated_proteins = SEQKIT_CONCAT (
        ch_proteins_to_concat.map { rec -> record(id: rec.id, fasta_files: rec.all_training_proteins) }
    )

    ch_single_protein_db = ch_single_protein_db.map { rec -> rec + record(fasta: rec.all_training_proteins[0]) }

    ch_no_protein_db.map { rec -> warnNoProteinsFound(rec) }

    // ----------------------------------------------------------
    // CHECK HEADERS OR WHOLE PROTEIN DB AND FIX THEM WHEN NECESSARY
    // ----------------------------------------------------------

    ch_proteins = CHECK_PROTEIN_FASTA(
        ch_single_protein_db.mix( ch_concatenated_proteins ),
        min_prot_db_seq_length
    )

    ch_proteins = ch_proteins.map { rec -> record(id: rec.id, proteins_fasta: rec.fasta) }
    
    emit:
    ch_input.join(ch_proteins, by: 'id')
}

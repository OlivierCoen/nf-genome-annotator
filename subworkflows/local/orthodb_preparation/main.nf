nextflow.enable.types = true

include { ORTHODB_GET_AVAILABLE_CLADES                          } from '../../../modules/local/orthodb/get_available_clades'
include { ORTHODB_GET_CLADE                                     } from '../../../modules/local/orthodb/get_clade'
include { ORTHODB_MAKECLADEDB                                   } from '../../../modules/local/orthodb/make_clade_db'


def getOrthoDBHash(rec){
    return "${rec.orthodb_clade}_${rec.orthodb_excluded_clades.join("_")}_${rec.orthodb_excluded_species.join("_")}".hashCode()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    species: String
    orthodb_excluded_clades: List<String>
    orthodb_excluded_species: List<String>
}

workflow ORTHODB_PREPARATION {

    take:
    ch_input: Channel<Input>
  

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FETCH LIST OF ORTHODB CLADES
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_available_orthodb_clades = ORTHODB_GET_AVAILABLE_CLADES([])

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SEARCH A MATCHING CLADE BETWEEN THE LIST OF AVAILABLE CLADES AND THE LINEAGE OF PARENT ORGANISMS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ch_orthodb_clade = ORTHODB_GET_CLADE( 
        ch_input.map { rec -> rec.species }.unique(),
        ch_available_orthodb_clades
    )

    ch_input = ch_input.join(ch_orthodb_clade, by: 'species')

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DOWNLOADS THE CORRESPONDING ORTHODB CLADE
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_input = ch_input.map { rec -> rec + record(orthodb_hash: getOrthoDBHash(rec)) }

    // Printing message for each sample where no orthodb clade could be found
    ch_input
        .filter { rec -> rec.orthodb_clade == null }
        .map { rec -> 
            println "No OrthoDB clade could be found for species ${rec.species}. Skipping preparation of OrthoDB proteins for this species"
        }

    ch_orthodo_cladedb = ORTHODB_MAKECLADEDB( 
        ch_input.filter { rec -> rec.orthodb_clade != null }
    )

    ch_input = ch_input
                .join(
                    ch_orthodo_cladedb.map { rec -> rec + record(orthodb_hash: getOrthoDBHash(rec)) },
                    by: 'orthodb_hash'
                )
                .map { rec -> rec.subMap(rec.keySet() - ['orthodb_hash']) } // remove 'orthodb_hash' entry

    emit:
    ch_input

}

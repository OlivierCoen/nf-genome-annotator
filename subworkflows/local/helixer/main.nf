nextflow.enable.types = true

include { HELIXER_GET_LINEAGE as GET_LINEAGE       } from '../../../modules/local/helixer/get_lineage'
include { HELIXER_HELIXER     as RUN               } from '../../../modules/local/helixer/helixer'
include { HELIXER_FETCHMODEL  as FETCH_MODEL       } from '../../../modules/local/helixer/fetch_model'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    species: String
    fasta: Path
}


workflow HELIXER {

    take:
    ch_input: Channel<Input>

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SEARCH A MATCHING LINEAGE BETWEEN THE LIST OF AVAILABLE HELIXER LINEAGES AND THE LINEAGE OF PARENT ORGANISMS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_helixer_lineage = GET_LINEAGE( 
        ch_input.map { rec -> rec.species }.unique()
    )
    ch_input = ch_input.join(ch_helixer_lineage, by: 'species')

    // Printing message for each sample where no helixer lineage could be found
    ch_input
        .filter { rec -> rec.helixer_lineage == null }
        .map { rec -> 
            println "No Helixer lineage could be found for species ${rec.species}. Skipping structural annotation for this species."
        }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DOWNLOAD MODEL
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_helixer_model_to_download = ch_input
                                        .filter { rec -> rec.helixer_lineage != null }
                                        .map { rec -> rec.helixer_lineage }
                                        .unique()
                                        
    ch_helixer_model = FETCH_MODEL( ch_helixer_model_to_download )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // RUN HELIXER
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_helixer_input = ch_input
                        .join( 
                            ch_helixer_model.map { rec -> record(helixer_lineage: rec.lineage, helixer_models_path: rec.models_path) } ,
                            by: 'helixer_lineage'
                        )
                        .map { rec -> record(
                            id: rec.id,
                            fasta: rec.fasta,
                            lineage: rec.helixer_lineage,
                            models_path: rec.helixer_models_path
                        ) }

    ch_helixer_out = RUN( ch_helixer_input )

    ch_helixer_out = ch_helixer_out
            .filter { rec -> 
                def gff_lines = rec.gff.splitCsv(sep: '\t').findAll { !it[0].startsWith('#') }
                if ( gff_lines.size() == 0 ) {
                    log.warn("Helixer failed to predict any gene model for sample ${rec.id}. Skipping this sample.")
                }
                gff_lines.size() > 0
            }
            .map { rec -> record(id: rec.id, structural_annotation: rec.gff) }

    // samples for which helixer_lineage is null or helixer output is empty
    // will not be proceeded in the subsequent steps
    emit:
    annotated = ch_input.join( ch_helixer_out, by: 'id' )
    

}


nextflow.enable.types = true

include { HELIXER_HELIXER as RUN                  } from '../../../modules/local/helixer/helixer'
include { HELIXER_FETCHMODEL as FETCH_MODEL       } from '../../../modules/local/helixer/fetch_model'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    helixer_lineage: String
    fasta: Path
}


workflow HELIXER {

    take:
    ch_input: Channel<Input>

    main:

    // Printing message for each sample where no helixer lineage could be found
    ch_input
        .filter { rec -> rec.helixer_lineage == null }
        .map { rec -> 
            println "No Helixer lineage could be found for sample ${rec.id}. Skipping structural annotation for this sample."
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

    ch_input = ch_input.join(
        ch_helixer_out.map { rec -> record(id: rec.id, structural_annotation: rec.gff) },
        by: 'id'
    )

    // samples for which helixer_lineage is null will not be proceeded in the subsequent steps
    emit:
    annotated = ch_input
    

}


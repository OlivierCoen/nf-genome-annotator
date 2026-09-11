nextflow.enable.types = true

include {                   } from '../../../modules/local/tiberius/tiberius'
include {        } from '../../../modules/local/tiberius/fetch_model'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    tiberius_lineage: String
    fasta: Path
}


workflow TIBERIUS {

    take:
    ch_input: Channel<Input>

    main:

    // Printing message for each sample where no tiberius lineage could be found
    ch_input
        .filter { rec -> rec.tiberius_lineage == null }
        .map { rec -> 
            println "No Tiberius lineage could be found for sample ${rec.id}. Skipping structural annotation for this sample."
        }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DOWNLOAD MODEL
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_helixer_model_to_download = ch_input
                                        .filter { rec -> rec.tiberius_lineage != null }
                                        .map { rec -> rec.tiberius_lineage }
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


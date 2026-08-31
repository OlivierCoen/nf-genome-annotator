nextflow.enable.types = true

include { HELIXER_HELIXER               } from '../../../modules/local/helixer/helixer'

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

    ch_helixer_input = ch_input
                        .filter { rec -> rec.helixer_lineage != null }
                        .map { rec -> record(
                            id: rec.id,
                            fasta: rec.fasta,
                            lineage: rec.helixer_lineage,
                            species: rec.species
                        ) }

    ch_helixer_out = HELIXER_HELIXER( ch_helixer_input )

    // samples for which helixer_lineage is null will not be proceeded in the subsequent steps
    emit:
    annotated = ch_input.join(ch_helixer_out, by: 'id')
    

}


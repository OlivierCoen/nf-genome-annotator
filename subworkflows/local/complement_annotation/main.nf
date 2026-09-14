nextflow.enable.types = true

include { AGAT_SPCOMPLEMENTANNOTATIONS       } from '../../../modules/local/agat/spcomplementannotations'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    structural_annotation: Path
    reference_gff: Path
}


workflow COMPLEMENT_ANNOTATION {

    take:
    ch_input: Channel<Input>
    complementation_use_ref_as_template: Boolean
    
    main:

    // ----------------------------------------------------------
    // PREPARE PROTEIN TRAINING SET FOR BRAKER
    // ----------------------------------------------------------

    ch_to_complement  = ch_input.filter{ rec -> rec.reference_gff != null }
    ch_leave_me_alone = ch_input.filter{ rec -> rec.reference_gff == null }

    // saving the uncomplemented annotation
    ch_to_complement = ch_to_complement.map { rec -> rec + record(uncomplemented_annotation: rec.structural_annotation) }

    ch_to_complement = ch_to_complement.map { rec ->
        if (complementation_use_ref_as_template) {
            record(
                id: rec.id, 
                ref_gff: rec.reference_gff, 
                other_gff: rec.structural_annotation
            )
        } else {
            record(
                id: rec.id, 
                ref_gff: rec.structural_annotation, 
                other_gff: rec.reference_gff
            ) 
        }
    }

    ch_complemented = AGAT_SPCOMPLEMENTANNOTATIONS( ch_to_complement )

    ch_complemented = ch_to_complement.join( 
        ch_complemented.map{ rec -> record(id: rec.id, structural_annotation: rec.gff)}, 
        by: 'id' 
    )


    emit:
    ch_complemented.mix( ch_leave_me_alone )

}

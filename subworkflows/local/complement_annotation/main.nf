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
    gff: Path
}


workflow COMPLEMENT_ANNOTATION {

    take:
    ch_input: Channel<Input>
    

    main:

    // ----------------------------------------------------------
    // PREPARE PROTEIN TRAINING SET FOR BRAKER
    // ----------------------------------------------------------

    ch_to_complement  = ch_input.filter{ rec -> rec.gff != null }
    ch_leave_me_alone = ch_input.filter{ rec -> rec.gff == null }

    // saving the uncomplemented annotation
    ch_to_complement = ch_to_complement.map { rec -> rec + record(uncomplemented_annotation: rec.structural_annotation) }

    ch_complemented = AGAT_SPCOMPLEMENTANNOTATIONS( 
        ch_to_complement.map{ rec -> record(
            id: rec.id, 
            ref_gff: rec.structural_annotation, 
            other_gff: rec.gff
        ) }
    )

    ch_complemented = ch_to_complement.join( 
        ch_complemented.map{ rec -> record(id: rec.id, annotation: rec.gff)}, 
        by: 'id' 
    )


    emit:
    complemented = ch_complemented.mix( ch_leave_me_alone )

}

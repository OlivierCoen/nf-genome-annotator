nextflow.enable.types = true

include { GFFREAD as FROM_MAIN_ANNOTATION             } from '../../../modules/local/gffread'
include { GFFREAD as FROM_OTHER_ANNOTATIONS           } from '../../../modules/local/gffread'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    gff: Path
    intermediate_annotations: Iterable<Path>
    alternative_annotations: Iterable<Path>
    previous_annotation: Path?
}

workflow EXTRACT_SEQUENCES {

    take:
    ch_input: Channel<Input>

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // EXTRACT PROTEOME FROM ACTUAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ch_proteome = FROM_MAIN_ANNOTATION( ch_input )
    
    ch_input = ch_input.join( 
        ch_proteome.map { rec -> record(id: rec.id, proteome: rec.proteins) }, 
        by: 'id'
    ) 

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GETTING PROTEOMES FOR ALL INTERMEDIATES + ALTERNATIVES [ + PREVIOUS ONE ] ) ANNOTATIONS
    // THIS WILL BE USED BY BUSCO
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ch_others_input = ch_input.flatMap { rec -> 
        def alt_gff = rec.alternative_annotations ?: []
        def interm_gffs = rec.intermediate_annotations ?: []
        def prev_gff = rec.previous_annotation ? [rec.previous_annotation]: []
        def other_gffs = alt_gff + interm_gffs + prev_gff
        other_gffs.collect { gff_file -> record(id: rec.id, gff: gff_file, fasta: rec.fasta) }
    }
                                 
    ch_other_proteomes = FROM_OTHER_ANNOTATIONS( ch_others_input )
    
    ch_other_proteomes = ch_other_proteomes.map { rec -> tuple(rec.id, rec.proteins) }
                            .groupTuple()
                            .map { id, fasta_list -> record(id: id, other_proteomes: fasta_list) }

    ch_input = ch_input.join( ch_other_proteomes, by: 'id' )

    emit:
    ch_input

}

nextflow.enable.types = true

include { AGAT_SPEXTRACTSEQUENCES as EXTRACT_PROTEOME             } from '../../../modules/local/agat/spextractsequences'
include { AGAT_SPEXTRACTSEQUENCES as EXTRACT_OTHER_PROTEOMES      } from '../../../modules/local/agat/spextractsequences'
include { AGAT_SPEXTRACTSEQUENCES as EXTRACT_CDS                  } from '../../../modules/local/agat/spextractsequences'

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
    codon_usage_id: Integer

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // EXTRACT PROTEOME FROM ACTUAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_proteome = EXTRACT_PROTEOME(
        ch_input,
        codon_usage_id
    )
    ch_input = ch_input
                .join( 
                    ch_proteome.map { rec -> record(id: rec.id, proteome: rec.extracted_fasta) }, 
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
        other_gffs.collect { gff -> record(id: rec.id, gff: gff, genome: rec.fasta) }
    }
                                    
    ch_other_proteomes = EXTRACT_OTHER_PROTEOMES( 
        ch_others_input,
        codon_usage_id
    ) 
    
    ch_other_proteomes = ch_other_proteomes.map { rec -> tuple(rec.id, rec.extracted.fasta) }
                            .groupTuple()
                            .map { id, fasta_list -> record(id: id, other_proteomes: fasta_list) }

    ch_input = ch_input.join( ch_other_proteomes, by: 'id' )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GETTING CDS SEQUENCES FOR THE MAIN ANNOTATION ONLY
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_cds = EXTRACT_CDS(
        ch_input.map { rec -> record(id: rec.id, gff: rec.gff, genome: rec.fasta) },
        null
    )
    ch_input = ch_input
                .join( 
                    ch_cds.map { rec -> record(id: rec.id, cds: rec.extracted_fasta) }, 
                    by: 'id'
                )

    emit:
    ch_input

}

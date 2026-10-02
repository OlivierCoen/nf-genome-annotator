nextflow.enable.types = true

include { MINIMAP2_INDEX                  } from '../../../modules/local/minimap2/index'
include { MINIMAP2_ALIGN                  } from '../../../modules/local/minimap2/align'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record MappingInput {
    id: String
    read_id: String
    reads: List<Path>
    fasta: Path
}

workflow FASTQ_ALIGN_MINIMAP2 {

    take:
    ch_input: Channel<MappingInput>

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // INDEX GENOME FOR MINIMAP2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_minimap2_index = MINIMAP2_INDEX(
        ch_input.map { rec -> rec.subMap(['id', 'fasta']) }.unique()
    )

    ch_input_with_index = ch_input.join(ch_minimap2_index, by: 'id')

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_aligned = MINIMAP2_ALIGN( 
        ch_input_with_index.join(ch_input_with_index, by: 'read_id')
    )
        
    emit:
    mapped = ch_input.join(ch_aligned, by: 'read_id')
}

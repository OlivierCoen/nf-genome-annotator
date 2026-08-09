nextflow.enable.types = true

include { HISAT2_EXTRACTSPLICESITES     } from '../../../modules/local/hisat2/extractsplicesites'
include { HISAT2_EXTRACTEXONS           } from '../../../modules/local/hisat2/extractexons'
include { HISAT2_BUILD                  } from '../../../modules/local/hisat2/build'
include { HISAT2_ALIGN                  } from '../../../modules/local/hisat2/align'

record MappingInput {
    id: String
    reads: List<Path>
    fasta: Path
    gtf: Path
}

workflow FASTQ_ALIGN_HISAT2 {

    take:
    ch_input: Channel<MappingInput>
    ignore_existing_gff_for_mapping: Boolean

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // EXTRACT SPLICE SITES AND EXONS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !ignore_existing_gff_for_mapping ) {

        HISAT2_EXTRACTSPLICESITES( ch_input )

        HISAT2_EXTRACTEXONS( ch_input )

        ch_input = ch_input
                    .join(HISAT2_EXTRACTSPLICESITES.out, by: 'id')
                    .join(HISAT2_EXTRACTEXONS.out, by: 'id')

    } else {
        ch_input = ch_input.map{ rec -> rec + record(splice_sites: null, exons: null) }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // INDEX GENOME FOR HISAT2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    HISAT2_BUILD( ch_input )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    HISAT2_ALIGN(
        ch_input.join(HISAT2_BUILD.out, by: 'id')
    )

    emit:
    mapped = ch_input.join(HISAT2_ALIGN, by: 'id')
}

nextflow.enable.types = true

include { HISAT2_EXTRACTSPLICESITES     } from '../../../modules/local/hisat2/extractsplicesites'
include { HISAT2_EXTRACTEXONS           } from '../../../modules/local/hisat2/extractexons'
include { HISAT2_BUILD                  } from '../../../modules/local/hisat2/build'
include { HISAT2_ALIGN                  } from '../../../modules/local/hisat2/align'

record MappingInput {
    id: String
    reads: List<Path>
    fasta: Path
    reference_gtf: Path
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

        ch_input_with_gtf = ch_input
                                .filter{ rec -> rec.reference_gtf != null }
                                .map{ rec -> record(id: rec.id, gtf: rec.reference_gtf) }
                                .unique()

        ch_splice_sites = HISAT2_EXTRACTSPLICESITES( ch_input_with_gtf )

        ch_exons = HISAT2_EXTRACTEXONS( ch_input_with_gtf )

        ch_input_with_gtf = ch_input_with_gtf
                                .join(ch_splice_sites, by: 'id')
                                .join(ch_exons, by: 'id')

        ch_input_without_gtf = ch_input
                                .filter{ rec -> rec.reference_gtf == null }
                                .map{ rec -> rec + record(splice_sites: null, exons: null) }

        ch_input = ch_input_with_gtf.mix( ch_input_without_gtf )

    } else {
        ch_input = ch_input.map{ rec -> rec + record(splice_sites: null, exons: null) }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // INDEX GENOME FOR HISAT2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_hisat2_index = HISAT2_BUILD( ch_input )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_aligned = HISAT2_ALIGN(
        ch_input.join(ch_hisat2_index, by: 'id')
    )

    emit:
    mapped = ch_input.join(ch_aligned, by: 'id')
}

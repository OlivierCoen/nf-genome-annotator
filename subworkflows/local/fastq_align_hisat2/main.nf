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

        ch_input_with_gtf = ch_input
                                .filter{ rec -> rec.gtf != null }
                                .map{ rec -> rec.subMap(['sample_id', 'gtf']) }
                                .unique()

        HISAT2_EXTRACTSPLICESITES( ch_input_with_gtf )

        HISAT2_EXTRACTEXONS( ch_input_with_gtf )

        ch_input_with_gtf = ch_input_with_gtf
                                .join(HISAT2_EXTRACTSPLICESITES.out, by: 'sample_id')
                                .join(HISAT2_EXTRACTEXONS.out, by: 'sample_id')

        ch_input_without_gtf = ch_input
                                .filter{ rec -> rec.gtf == null }
                                .map{ rec -> rec + record(splice_sites: null, exons: null) }

        ch_input = ch_input_with_gtf.mix( ch_input_without_gtf )

    } else {
        ch_input = ch_input.map{ rec -> rec + record(splice_sites: null, exons: null) }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // INDEX GENOME FOR HISAT2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_fasta_to_build = ch_input
                            .map { rec -> rec.subMap(['sample_id', 'fasta', 'splice_sites', 'exons']) }
                            .unique()

    HISAT2_BUILD( ch_fasta_to_build )

    ch_input = ch_input.join(HISAT2_BUILD.out, by: 'sample_id')

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    HISAT2_ALIGN(
        ch_input.join(HISAT2_BUILD.out, by: 'sample_id')
    )

    emit:
    mapped = ch_input.join(HISAT2_ALIGN.out, by: 'id')
}

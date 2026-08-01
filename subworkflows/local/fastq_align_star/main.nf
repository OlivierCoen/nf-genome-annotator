nextflow.enable.types = true

include { STAR_GENOMEGENERATE            } from '../../../modules/local/star/genomegenerate'
include { STAR_ALIGN                     } from '../../../modules/local/star/align'

record MappingInput {
    id: String
    reads: List<Path>
    fasta: Path
    gtf: Path
}

workflow FASTQ_ALIGN_STAR {

    take:
    ch_input: Channel<MappingInput>
    ignore_existing_gff_for_mapping: Boolean

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // INDEX GENOME FOR STAR
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    STAR_GENOMEGENERATE(
        ch_input,
        ignore_existing_gff_for_mapping
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    STAR_ALIGN(
        ch_input.join(STAR_GENOMEGENERATE.out, by: 'id'),
        ignore_existing_gff_for_mapping
    )

    emit:
    mapped = ch_input.join(STAR_ALIGN.out, by: 'id')

}

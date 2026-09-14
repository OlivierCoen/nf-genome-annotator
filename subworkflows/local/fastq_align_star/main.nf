nextflow.enable.types = true

include { STAR_GENOMEGENERATE            } from '../../../modules/local/star/genomegenerate'
include { STAR_ALIGN                     } from '../../../modules/local/star/align'

record MappingInput {
    id: String
    reads: List<Path>
    fasta: Path
    reference_gtf: Path
}

workflow FASTQ_ALIGN_STAR {

    take:
    ch_input: Channel<MappingInput>
    ignore_existing_gff_for_mapping: Boolean

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // INDEX GENOME FOR STAR
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_star_index = STAR_GENOMEGENERATE(
        ch_input.map { rec -> record(id: rec.id, fasta: rec.fasta, gtf: rec.reference_gtf) },
        ignore_existing_gff_for_mapping
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_aligned = STAR_ALIGN(
        ch_input.join(ch_star_index, by: 'id')
    )

    emit:
    mapped = ch_input.join(ch_aligned, by: 'id')

}

include { STAR_GENOMEGENERATE            } from '../../../modules/local/star/genomegenerate'
include { STAR_ALIGN                     } from '../../../modules/local/star/align'



workflow FASTQ_ALIGN_STAR {

    take:
    ch_genome
    ch_reads
    ch_gtf
    star_ignore_existing_gtf

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // INDEX GENOME FOR STAR
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_star_genomegenerate_input = ch_genome
                                    .join( ch_gtf, remainder: true ) // gives objects even if no corresponding gtf
                                    .map{
                                        meta, genome, gtf ->
                                            [ meta, genome, genome ]
                                    }
                                    .view{ v -> "genome $v"}

    STAR_GENOMEGENERATE(
        ch_star_genomegenerate_input,
        star_ignore_existing_gtf
    )
    ch_index = STAR_GENOMEGENERATE.out.index

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (star_ignore_existing_gtf) {
        ch_star_input = ch_reads
                            .cross( ch_index ) { v -> v[0].id } // match only on id, ignore single_end
                            .map{ meta, reads, index -> [ meta, reads, index, [] ] }
    } else {
        ch_star_input = ch_reads
                            .cross( ch_index ) { v -> v[0].id }
                            .cross( ch_gtf ) { v -> v[0].id }
    }

    STAR_ALIGN(
        ch_star_input,
        star_ignore_existing_gtf
    )

    emit:
    bam                 = STAR_ALIGN.out.bam // channel: [ val(meta), path(bam) ]

}

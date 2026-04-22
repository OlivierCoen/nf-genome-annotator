include { FASTP                          } from '../../../modules/nf-core/fastp'
include { STAR_GENOMEGENERATE            } from '../../../modules/local/star/genomegenerate'
include { STAR_ALIGN                     } from '../../../modules/local/star/align'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow MAP_RNASEQ_READS {

    take:
    ch_genome
    ch_reads
    ch_gtf
    star_ignore_existing_gtf

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PREPROCESS READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    FASTP(
        ch_reads.map{ meta, reads -> [ meta, reads, [] ] },
        [], [], []
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // INDEX GENOME FOR STAR
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_star_genomegenerate_input = star_ignore_existing_gtf ?
                                    ch_genome.map{ meta, genome -> [ meta, genome, [] ] } :
                                    ch_genome.join( ch_gtf )

    STAR_GENOMEGENERATE( ch_star_genomegenerate_input )
    ch_index = STAR_GENOMEGENERATE.out.index

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_star_input = star_ignore_existing_gtf ?
                        ch_reads.join( ch_index ).map{ meta, reads, index -> [ meta, reads, index, [] ] } :
                        ch_reads.join( ch_index ).join( ch_gtf )

    STAR_ALIGN( ch_star_input )

    // keeping only meta id in meta map
    ch_bam = STAR_ALIGN.out.bam
                .map { meta, bam -> [ [ id: meta.id ], bam ] }

    emit:
    bam = ch_bam

}

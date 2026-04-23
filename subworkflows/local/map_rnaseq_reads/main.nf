include { FASTP                          } from '../../../modules/nf-core/fastp'
include { STAR_GENOMEGENERATE            } from '../../../modules/local/star/genomegenerate'
include { STAR_ALIGN                     } from '../../../modules/local/star/align'

include { FASTQ_FASTQC_UMITOOLS_FASTP    } from '../../nf-core/fastq_fastqc_umitools_fastp'



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
    skip_fastqc
    skip_umi_extract
    skip_trimming
    star_ignore_existing_gtf

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FASTQC & FASTP
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    FASTQ_FASTQC_UMITOOLS_FASTP(
        ch_reads.map{ meta, reads -> [ meta, reads, [] ] }, // add empty spot for adapter
        skip_fastqc,
        with_umi=true,
        skip_umi_extract,
        umi_discard_read=0,
        skip_trimming,
        save_trimmed_fail=false,
        save_merged=false,
        min_trimmed_reads=0

    )

    ch_reads = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // INDEX GENOME FOR STAR
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_star_genomegenerate_input = ch_genome
                                    .join( ch_gtf, remainder: true ) // gives objects even if no corresponding gtf
                                    .map{
                                        meta, genome, gtf ->
                                            if ( star_ignore_existing_gtf ) {
                                                [ meta, genome, [] ]
                                            } else {
                                                [ meta, genome, gtf ]
                                            }
                                    }
                                    .view()

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

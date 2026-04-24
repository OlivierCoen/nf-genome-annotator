include { FASTQ_FASTQC_UMITOOLS_FASTP    } from '../../nf-core/fastq_fastqc_umitools_fastp'
include { FASTQ_ALIGN_STAR               } from '../fastq_align_star'
include { BAM_SORT_INDEX_STATS           } from '../bam_sort_index_stats'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow MAP_TO_GENOME_SORT_INDEX {

    take:
    ch_genome
    ch_reads
    ch_bam
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
    // MAP READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    FASTQ_ALIGN_STAR(
        ch_genome,
        ch_reads,
        ch_gtf,
        star_ignore_existing_gtf
    )

    ch_aligned_bam = FASTQ_ALIGN_STAR.out.bam
                        .map { meta, bam -> [ [ id: meta.id ], bam ] } // keeping only meta.id in meta map

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAPPING STATISTICS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    BAM_SORT_INDEX_STATS(
        ch_bam.mix( ch_aligned_bam ),
        ch_genome
    )


    emit:
    bam_bai = BAM_SORT_INDEX_STATS.out.bam_bai

}

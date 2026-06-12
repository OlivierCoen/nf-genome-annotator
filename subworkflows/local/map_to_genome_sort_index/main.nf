include { FASTQ_FASTQC_UMITOOLS_FASTP    } from '../../nf-core/fastq_fastqc_umitools_fastp'
include { FASTQ_ALIGN_HISAT2             } from '../fastq_align_hisat2'
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
    rnaseq_mapper
    ignore_existing_gtf_for_mapping

    main:

    // get only genomes that need to be built (genomes for which there are reads)
    ch_genome_for_mapping = ch_reads
                            .cross( ch_genome ) { v -> v[0][0] } // match only on id
                            .map{ // [[meta, reads], [meta2, index]]
                                part1, genome_part -> 
                                    def meta = genome_part[0]
                                    def genome = genome_part[1]
                                    [ meta, genome ]
                            }

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

    if (rnaseq_mapper == 'hisat2') {

        FASTQ_ALIGN_HISAT2(
            ch_genome_for_mapping,
            ch_reads,
            ch_gtf,
            ignore_existing_gtf_for_mapping
        )
    
        ch_aligned_bam = FASTQ_ALIGN_HISAT2.out.bam
    
    } else if (rnaseq_mapper == 'star') {

        FASTQ_ALIGN_STAR(
            ch_genome_for_mapping,
            ch_reads,
            ch_gtf,
            ignore_existing_gtf_for_mapping
        )

        ch_aligned_bam = FASTQ_ALIGN_STAR.out.bam

    }

    ch_aligned_bam = ch_aligned_bam
                        .map { meta, bam -> [ [ id: meta.id ], bam ] } // keeping only meta.id in meta map

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAPPING STATISTICS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    BAM_SORT_INDEX_STATS(
        ch_bam.mix( ch_aligned_bam ),
        ch_genome_for_mapping
    )


    emit:
    bam_bai = BAM_SORT_INDEX_STATS.out.bam_bai

}

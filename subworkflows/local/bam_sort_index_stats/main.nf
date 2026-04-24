include { SAMTOOLS_FAIDX                           } from '../../../modules/local/samtools/faidx'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_INDEX     } from '../../../modules/local/samtools/sort'
include { SAMTOOLS_STATS                           } from '../../../modules/local/samtools/stats'
include { SAMTOOLS_IDXSTATS                        } from '../../../modules/nf-core/samtools/idxstats'
include { SAMTOOLS_FLAGSTAT                        } from '../../../modules/nf-core/samtools/flagstat'

workflow BAM_SORT_INDEX_STATS {
    take:
    ch_bam // channel: [ val(meta), path(bam) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]

    main:

    SAMTOOLS_FAIDX( ch_fasta )
    ch_fasta_fai = ch_fasta.join( SAMTOOLS_FAIDX.out.fai )

    SAMTOOLS_SORT_INDEX( ch_bam )
    ch_bam_bai = SAMTOOLS_SORT_INDEX.out.bam
                    .join( SAMTOOLS_SORT_INDEX.out.index )

    SAMTOOLS_STATS(
        ch_bam_bai.join( ch_fasta_fai )
    )

    SAMTOOLS_FLAGSTAT(ch_bam_bai)

    SAMTOOLS_IDXSTATS(ch_bam_bai)

    emit:
    bam_bai  = ch_bam_bai
    stats    = SAMTOOLS_STATS.out.stats // channel: [ val(meta), path(stats) ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), path(idxstats) ]
}

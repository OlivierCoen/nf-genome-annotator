nextflow.enable.types = true

include { SAMTOOLS_FAIDX                           } from '../../../modules/local/samtools/faidx'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_INDEX     } from '../../../modules/local/samtools/sort'
include { SAMTOOLS_STATS                           } from '../../../modules/local/samtools/stats'
include { SAMTOOLS_IDXSTATS                        } from '../../../modules/local/samtools/idxstats'
include { SAMTOOLS_FLAGSTAT                        } from '../../../modules/local/samtools/flagstat'

record Bams {
    id: String
    bam: Iterable<Path>
}

workflow BAM_SORT_INDEX_STATS {

    take:
    ch_input: Channel<Bams>

    main:

    // ------------------------------------------------------------------------------------
    // INDEX FASTA
    // ------------------------------------------------------------------------------------

    SAMTOOLS_FAIDX(
        ch_input.map { rec -> rec.subMap(['id', 'fasta']) }.unique()
    )
    ch_input = ch_input.join( SAMTOOLS_FAIDX.out, by: 'id' )

    // ------------------------------------------------------------------------------------
    // SORT BAMS AND MAKE INDEX
    // ------------------------------------------------------------------------------------

    ch_bam = ch_input.flatMap { rec ->
        rec.bams.collect{ bam ->
            record(
                sample_id: rec.id,
                id: bam.baseName,
                bam: bam,
                fasta: rec.fasta,
                fai: rec.fai
            ) }
    }

    SAMTOOLS_SORT_INDEX( ch_bam )

    ch_bam = ch_bam.join( SAMTOOLS_SORT_INDEX.out, by: 'id' )

    // ------------------------------------------------------------------------------------
    // MAPPING STATS
    // ------------------------------------------------------------------------------------

    SAMTOOLS_STATS( ch_bam )

    SAMTOOLS_FLAGSTAT( ch_bam )

    SAMTOOLS_IDXSTATS( ch_bam )

    // ------------------------------------------------------------------------------------
    // ASSOCIATE SORTED BAM TO ORIGINAL DATA
    // ------------------------------------------------------------------------------------

    ch_bams = ch_bam
                .map { rec -> tuple( rec.sample_id, record(bam: rec.bam, bai: rec.bai) ) }
                .groupTuple()
                .map { id, rec_list -> record(id: id, mappings: rec_list) }

    emit:
    sorted_indexed = ch_input.join( ch_bams, by: 'id' )

}

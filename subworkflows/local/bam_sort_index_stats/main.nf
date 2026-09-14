nextflow.enable.types = true

include { SAMTOOLS_FAIDX                           } from '../../../modules/local/samtools/faidx'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_INDEX     } from '../../../modules/local/samtools/sort'
include { SAMTOOLS_STATS                           } from '../../../modules/local/samtools/stats'
include { SAMTOOLS_IDXSTATS                        } from '../../../modules/local/samtools/idxstats'
include { SAMTOOLS_FLAGSTAT                        } from '../../../modules/local/samtools/flagstat'

record Bams {
    id: String
    short_read_bams: Set<Path>
}

workflow BAM_SORT_INDEX_STATS {

    take:
    ch_input: Channel<Bams>

    main:

    // ------------------------------------------------------------------------------------
    // INDEX FASTA
    // ------------------------------------------------------------------------------------

    ch_fai = SAMTOOLS_FAIDX(
        ch_input.map { rec -> rec.subMap(['id', 'fasta']) }.unique()
    )
    
    ch_input = ch_input.join( ch_fai, by: 'id' )

    // ------------------------------------------------------------------------------------
    // SORT BAMS AND MAKE INDEX
    // ------------------------------------------------------------------------------------

    ch_bam = ch_input.flatMap { rec ->
        rec.short_read_bams.collect{ bam ->
            record(
                id: rec.id,
                bam: bam,
                fasta: rec.fasta,
                fai: rec.fai
            ) }
    }

    ch_sorted_bam_bai = SAMTOOLS_SORT_INDEX( ch_bam )

    ch_bam_bai = ch_bam.join( ch_sorted_bam_bai, by: 'id' )

    // ------------------------------------------------------------------------------------
    // MAPPING STATS
    // ------------------------------------------------------------------------------------

    SAMTOOLS_STATS( ch_bam_bai )

    SAMTOOLS_FLAGSTAT( ch_bam_bai )

    SAMTOOLS_IDXSTATS( ch_bam_bai )

    // ------------------------------------------------------------------------------------
    // ASSOCIATE SORTED BAM TO ORIGINAL DATA
    // ------------------------------------------------------------------------------------

    ch_mappings = ch_bam_bai
                    .map { rec -> tuple( rec.id, record(bam: rec.bam, bai: rec.bai) ) }
                    .groupTuple()
                    .map { id, bam_bai_list -> record(
                        id: id, 
                        short_read_sorted_bams_bais: bam_bai_list
                    ) }

    emit:
    ch_input.join( ch_mappings, by: 'id' )

}

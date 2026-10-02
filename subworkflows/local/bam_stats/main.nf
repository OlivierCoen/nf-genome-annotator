nextflow.enable.types = true

include { SAMTOOLS_INDEX                           } from '../../../modules/local/samtools/index'
include { SAMTOOLS_STATS                           } from '../../../modules/local/samtools/stats'
include { SAMTOOLS_IDXSTATS                        } from '../../../modules/local/samtools/idxstats'
include { SAMTOOLS_FLAGSTAT                        } from '../../../modules/local/samtools/flagstat'

record Mapping {
    type: String
    bam: Path
}

record Input {
    id: String
    mappings: Set<Mapping>
}

workflow BAM_STATS {

    take:
    ch_input: Channel<Input>

    main:

    // ------------------------------------------------------------------------------------
    // SORT BAMS AND MAKE INDEX
    // ------------------------------------------------------------------------------------

    ch_bam = ch_input.flatMap { rec ->
        rec.bams.collect{ bam -> record(id: rec.id, bam: bam) }
    }

    ch_index_out = SAMTOOLS_INDEX( ch_bam )

    ch_bam_bai = ch_bam.join( ch_index_out, by: 'id' )

    // ------------------------------------------------------------------------------------
    // MAPPING STATS
    // ------------------------------------------------------------------------------------

    SAMTOOLS_STATS( ch_bam_bai )

    SAMTOOLS_FLAGSTAT( ch_bam_bai )

    SAMTOOLS_IDXSTATS( ch_bam_bai )

}

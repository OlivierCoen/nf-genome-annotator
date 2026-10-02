nextflow.enable.types = true

include { SAMTOOLS_FAIDX                           } from '../../../modules/local/samtools/faidx'
include { SAMTOOLS_CHECK_SORTED_BY_COORDINATE      } from '../../../modules/local/samtools/check_sorted_by_coordinate'
include { SAMTOOLS_INDEX                           } from '../../../modules/local/samtools/index'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_INDEX     } from '../../../modules/local/samtools/sort'
include { SAMTOOLS_STATS                           } from '../../../modules/local/samtools/stats'
include { SAMTOOLS_IDXSTATS                        } from '../../../modules/local/samtools/idxstats'
include { SAMTOOLS_FLAGSTAT                        } from '../../../modules/local/samtools/flagstat'

record Mapping {
    type: String
    bam: Path
}

record Input {
    id: String
    fasta: Path
    mappings: Set<Mapping>
}

workflow BAM_SORT_INDEX_STATS {

    take:
    ch_input: Channel<Input>

    main:

    // ------------------------------------------------------------------------------------
    // INDEX FASTA
    // ------------------------------------------------------------------------------------

    ch_fai = SAMTOOLS_FAIDX(
        ch_input.map { rec -> rec.subMap(['id', 'fasta']) }.unique()
    )
    ch_input.view()
    ch_input = ch_input.join( ch_fai, by: 'id' )

    // ------------------------------------------------------------------------------------
    // CHECK IF BAMS ARE SORTED
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

    ch_sorted_status = SAMTOOLS_CHECK_SORTED_BY_COORDINATE( ch_bam )
    ch_bam_with_status = ch_bam.join(ch_sorted_status, by: 'read_id')

    // ------------------------------------------------------------------------------------
    // SORT BAMS IF NEEDED + MAKE INDEX
    // ------------------------------------------------------------------------------------


    ch_sorted_bam   = ch_bam_with_status.filter { rec -> rec.sorted == true }
    ch_unsorted_bam = ch_bam_with_status.filter { rec -> rec.sorted == false }

    ch_index_out      = SAMTOOLS_INDEX( ch_sorted_bam )
    ch_sort_index_out = SAMTOOLS_SORT_INDEX( ch_unsorted_bam )
    
    ch_sorted_bam_bai = ch_index_out.mix( ch_sort_index_out )

    ch_bam.view{ v-> "before $v"}
    ch_bam_bai = ch_bam.join( ch_sorted_bam_bai, by: 'id' )
    ch_bam_bai.view{ v-> "after $v"}

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
                        short_read_mappings: bam_bai_list
                    ) }

    emit:
    ch_input.join( ch_mappings, by: 'id' )

}

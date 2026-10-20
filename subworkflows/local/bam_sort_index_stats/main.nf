nextflow.enable.types = true

include { SAMTOOLS_FAIDX                           } from '../../../modules/local/samtools/faidx'
include { SAMTOOLS_INDEX                           } from '../../../modules/local/samtools/index'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_INDEX     } from '../../../modules/local/samtools/sort'
include { SAMTOOLS_STATS                           } from '../../../modules/local/samtools/stats'
include { SAMTOOLS_IDXSTATS                        } from '../../../modules/local/samtools/idxstats'
include { SAMTOOLS_FLAGSTAT                        } from '../../../modules/local/samtools/flagstat'

record Bam {
    id: String
    type: String
    bam: Path
}

record Input {
    id: String
    fasta: Path
    bams: Set<Bam>
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
    ch_input = ch_input.join( ch_fai, by: 'id' )

    // ------------------------------------------------------------------------------------
    // CHECK IF BAMS ARE SORTED
    // ------------------------------------------------------------------------------------
   
    ch_bam = ch_input.flatMap { rec ->
        rec.bams.collect{ subrec ->
            record(
                id: rec.id,
                read_id: subrec.read_id,
                bam: subrec.bam,
                bai: subrec.bai,
                type: subrec.type,
                fasta: rec.fasta,
                fai: rec.fai,
                user_supplied: subrec.user_supplied
            ) }
    }

    // ------------------------------------------------------------------------------------
    // SORT BAMS IF NEEDED + MAKE INDEX
    // AT THIS POINT, ONLY USER SUPPLIED BAMS SHOULD NEED SORTING
    // ------------------------------------------------------------------------------------

    ch_user_supplied_bam = ch_bam.filter { rec -> rec.user_supplied == true }
    ch_new_bam           = ch_bam.filter { rec -> rec.user_supplied != true }
    
    ch_index_out      = SAMTOOLS_INDEX( ch_new_bam )
    ch_sort_index_out = SAMTOOLS_SORT_INDEX( ch_user_supplied_bam )
    
    ch_sorted_bam_bai = ch_index_out.mix( ch_sort_index_out )

    ch_bam_bai = ch_bam.join( ch_sorted_bam_bai, by: 'read_id' )

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
                    .map { rec -> tuple( rec.id, record(bam: rec.bam, bai: rec.bai, type: rec.type) ) }
                    .groupBy()
                    .map { id, records -> record(id: id, mappings: records) }

    emit:
    ch_input.join( ch_mappings, by: 'id' )

}

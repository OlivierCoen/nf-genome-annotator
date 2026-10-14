nextflow.enable.types = true

include { SAMTOOLS_MERGE as MERGE_SHORT_READ_BAMS               } from '../../../modules/local/samtools/merge'
include { SAMTOOLS_MERGE as MERGE_LONG_READ_BAMS                } from '../../../modules/local/samtools/merge'
include { BRAKER3                                               } from '../../../modules/local/braker3'
include { BRAKER3 as BRAKER3_ISOSEQ                             } from '../../../modules/local/braker3'
include { TSEBRA_TSEBRA as TSEBRA                               } from '../../../modules/local/tsebra/tsebra'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Mapping {
    bam: Path
    bai: Path
}

record Input {
    id: String
    fasta: Path
    species: String
    short_read_mappings: List<Mapping>
    long_read_mappings: List<Mapping>
    proteins_fasta: Path
    tsebra_gtfs: Set<Path>
    tsebra_hintsfiles: Set<Path>
}


workflow BRAKER3_PIPELINE {

    take:
    ch_input: Channel<Input>

    main:

    ch_has_mappings   = ch_input.filter{ rec -> rec.short_read_mappings != null }
    ch_leave_me_alone = ch_input.filter{ rec -> rec.short_read_mappings == null }
    
    ch_merge_me       = ch_has_mappings.filter{ rec -> rec.short_read_mappings.size() > 1 }
    ch_not_no_merge   = ch_has_mappings.filter{ rec -> rec.short_read_mappings.size() <= 1 }

    // ----------------------------------------------------------
    // MERGE MULTIPLE BAM FILES INTO A SINGLE BAM WHEN NECESSARY
    // ----------------------------------------------------------

    ch_samtools_merge_input = ch_merge_me.map{ rec ->
        def bams = rec.short_read_mappings.collect { r -> r.bam }
        def bais = rec.short_read_mappings.collect { r -> r.bai }
        record(id: rec.id, bams: bams, bais: bais)
    }

    ch_merged = MERGE_SHORT_READ_BAMS( ch_samtools_merge_input )

    ch_merge_me       = ch_merge_me.join( ch_merged, by: 'id' )
    ch_not_no_merge   = ch_not_no_merge.map { rec -> rec + record(bam: rec.short_read_mappings[0].bam) }

    ch_braker_input = ch_leave_me_alone.mix( ch_not_no_merge ).mix( ch_merge_me )

    // ----------------------------------------------------------
    // RUN BRAKER3
    // ----------------------------------------------------------

    ch_braker_out = BRAKER3( ch_braker_input )

    ch_braker_out = ch_braker_out
            .filter { rec -> 
                if ( rec.braker_gtf == null ) {
                    log.warn("Braker3 failed to predict any gene model for sample ${rec.id}. Skipping this sample.")
                }
                rec.braker_gtf != null
            }

    ch_input = ch_input.join(ch_braker_out, by: 'id')

    // ----------------------------------------------------------
    // MERGE ANNOTATIONS WHEN NECESSARY
    // ----------------------------------------------------------

    // separate inputs that need to be merged from the rest
    // normally, tsebra_gtfs and tsebra_hintsfiles should be already both present or both absent
    ch_to_merge_with_tsebra = ch_input.filter{ rec -> rec.tsebra_gtfs.size() > 0  && rec.tsebra_hintsfiles.size() > 0 }
    ch_not_to_merge         = ch_input.filter{ rec -> rec.tsebra_gtfs.size() == 0 || rec.tsebra_hintsfiles.size() == 0 }

    ch_tsebra_input = ch_to_merge_with_tsebra.map{ rec ->
        record(
            id: rec.id,
            gtfs: [rec.braker_gtf] + rec.tsebra_gtfs,
            hintsfiles: [rec.braker_hintsfile] + rec.tsebra_hintsfiles
        )
    }

    ch_tsebra_out = TSEBRA( ch_tsebra_input )

    ch_merged = ch_to_merge_with_tsebra.join(ch_tsebra_out, by: 'id')

    // ----------------------------------------------------------
    // MIXING MERGED AND NOT MERGED
    // ----------------------------------------------------------

    ch_merged       = ch_merged.map      { rec -> rec + record(structural_annotation: rec.merged_gtf) }
    ch_not_to_merge = ch_not_to_merge.map{ rec -> rec + record(structural_annotation: rec.braker_gtf) }

    emit:
    ch_merged.mix( ch_not_to_merge )

}

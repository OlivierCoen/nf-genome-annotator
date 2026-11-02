nextflow.enable.types = true

include { SAMTOOLS_MERGE as MERGE_BAMS                          } from '../../../modules/local/samtools/merge'
include { SAMTOOLS_MERGE as MERGE_LONG_READ_BAMS                } from '../../../modules/local/samtools/merge'
include { BRAKER3 as BRAKER                                     } from '../../../modules/local/braker3'
include { BRAKER3 as BRAKER_ISOSEQ                              } from '../../../modules/local/braker3'
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
    mappings: List<Mapping>
    long_read_mappings: List<Mapping>
    proteins_fasta: Path
    tsebra_gtfs: Set<Path>
    tsebra_hintsfiles: Set<Path>
}


workflow BRAKER3 {

    take:
    ch_input: Channel<Input>

    main:

    // ----------------------------------------------------------
    // SHORT READ BAMS
    // ----------------------------------------------------------

    // ----------------------------------------------------------
    // MERGE MULTIPLE SHORT READ BAM FILES INTO A SINGLE BAM WHEN NECESSARY
    // ----------------------------------------------------------
ch_input.view()
    ch_has_mappings   = ch_input.filter{ rec -> rec.short_read_mappings != null }
    ch_leave_me_alone = ch_input.filter{ rec -> rec.short_read_mappings == null }
    
    ch_merge_me       = ch_has_mappings.filter{ rec -> rec.short_read_mappings.size() > 1 }
    ch_not_no_merge   = ch_has_mappings.filter{ rec -> rec.short_read_mappings.size() <= 1 }

    ch_samtools_merge_input = ch_merge_me.map{ rec ->
        def bams = rec.short_read_mappings.collect { r -> r.bam }
        def bais = rec.short_read_mappings.collect { r -> r.bai }
        record(id: rec.id, bams: bams, bais: bais)
    }

    ch_merged = MERGE_BAMS( ch_samtools_merge_input )

    ch_merge_me      = ch_merge_me.join( ch_merged, by: 'id' )
    ch_not_no_merge  = ch_not_no_merge.map { rec -> rec + record(bam: rec.short_read_mappings[0].bam) }

    ch_braker_input = ch_leave_me_alone
                        .mix( ch_not_no_merge )
                        .mix( ch_merge_me )

    // ----------------------------------------------------------
    // RUN BRAKER3
    // ----------------------------------------------------------

    ch_braker_out = BRAKER( ch_braker_input )

    ch_braker_out = ch_braker_out
            .filter { rec -> 
                if ( rec.braker_gtf == null ) {
                    log.warn("Braker3 failed to predict any gene model for sample ${rec.id}. Skipping this sample.")
                }
                rec.braker_gtf != null
            }

    
    // ----------------------------------------------------------
    // LONG READ BAMS
    // ----------------------------------------------------------

    // ----------------------------------------------------------
    // MERGE MULTIPLE BAM FILES INTO A SINGLE BAM WHEN NECESSARY
    // ----------------------------------------------------------

    ch_has_isoseq_mappings   = ch_input.filter{ rec -> rec.long_read_mappings != null }
    
    ch_isoseq_merge_me       = ch_has_isoseq_mappings.filter{ rec -> rec.long_read_mappings.size() > 1 }
    ch_isoseq_not_no_merge   = ch_has_isoseq_mappings.filter{ rec -> rec.long_read_mappings.size() <= 1 }

    ch_isoseq_samtools_merge_input = ch_isoseq_merge_me.map{ rec ->
        def bams = rec.long_read_mappings.collect { r -> r.bam }
        def bais = rec.long_read_mappings.collect { r -> r.bai }
        record(id: rec.id, bams: bams, bais: bais)
    }

    ch_isoseq_merged = MERGE_LONG_READ_BAMS( ch_isoseq_samtools_merge_input )

    ch_isoseq_merge_me     = ch_isoseq_merge_me.join( ch_isoseq_merged, by: 'id' )
    ch_isoseq_not_no_merge = ch_isoseq_not_no_merge.map { rec -> rec + record(bam: rec.long_read_mappings[0].bam) }

    ch_isoseq_braker_input = ch_isoseq_not_no_merge.mix( ch_isoseq_merge_me )

    // ----------------------------------------------------------
    // RUN BRAKER3
    // ----------------------------------------------------------

    ch_isoseq_braker_out = BRAKER_ISOSEQ( ch_isoseq_braker_input )

    ch_isoseq_braker_out = ch_isoseq_braker_out
            .filter { rec -> 
                if ( rec.braker_gtf == null ) {
                    log.warn("Braker3 Isoseq failed to predict any gene model for sample ${rec.id}. Skipping this sample.")
                }
                rec.braker_gtf != null
            }
            .map { rec -> record(
                    id: rec.id,
                    isoseq_braker_gtf: rec.braker_gtf,
                    isoseq_braker_hintsfile: rec.braker_hintsfile
                )
            }



    // ----------------------------------------------------------
    // ASSOCIATING BRAKER RESULTS TO INPUT
    // ----------------------------------------------------------

    ch_input = ch_input
                    .join(ch_braker_out, by: 'id', remainder: true)
                    .join(ch_isoseq_braker_out, by: 'id', remainder: true)
                    .filter{ rec -> rec.braker_gtf != null || rec.isoseq_braker_gtf != null } // keeping only ercords for which at least a braker gtf is present
                    .map { rec ->
                        def newrec = rec.subMap(rec.keySet() - ['braker_gtf', 'braker_hintsfile', 'isoseq_braker_gtf', 'isoseq_braker_hintsfile']) + record(
                            tsebra_gtfs: rec.supplied_tsebra_gtfs + [rec.braker_gtf, rec.isoseq_braker_gtf].findAll(),
                            tsebra_hintsfiles: rec.supplied_tsebra_hintsfiles + [rec.braker_hintsfile, rec.isoseq_braker_hintsfile].findAll()
                        )
                        // checking that GTFs and Hintfiles for TSEBRA are the same number
                        if (newrec.tsebra_gtfs.size() != newrec.tsebra_hintsfiles.size() ) {
                            error("Different number of GTFs and Hintfiles given to TSEBRA. Please check the output of Braker runs.")
                        }
                        newrec
                    }
                   
    // ----------------------------------------------------------
    // MERGE ANNOTATIONS WHEN NECESSARY
    // ----------------------------------------------------------

    // separate inputs that need to be merged from the rest
    // normally, tsebra_gtfs and tsebra_hintsfiles should be already both present or both absent
    ch_not_to_merge         = ch_input.filter{ rec -> rec.tsebra_gtfs.size() == 0 }
    ch_to_merge_with_tsebra = ch_input.filter{ rec -> rec.tsebra_gtfs.size() > 0  }
                                
    ch_tsebra_input = ch_to_merge_with_tsebra.map{ rec ->
        record(
            id: rec.id,
            gtfs: rec.tsebra_gtfs,
            hintsfiles: rec.tsebra_hintsfiles
        )
    }

    ch_tsebra_out = TSEBRA( ch_tsebra_input )

    ch_merged = ch_to_merge_with_tsebra.join(ch_tsebra_out, by: 'id')

    // ----------------------------------------------------------
    // MIXING MERGED AND NOT MERGED
    // ----------------------------------------------------------

    ch_merged       = ch_merged.map      { rec -> rec + record(structural_annotation: rec.merged_gtf) }
    // for records that did not need merging, the only gtf in tsebra_gtfs is the one issued by the only Braker run that was launched (or that worked)
    ch_not_to_merge = ch_not_to_merge.map{ rec -> rec + record(structural_annotation: rec.tsebra_gtfs[0]) }

    emit:
    ch_merged.mix( ch_not_to_merge )

}

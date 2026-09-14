nextflow.enable.types = true

include { TRAINING_PROTEIN_PREPARATION                          } from '../training_protein_preparation'
include { BRAKER3                                               } from '../../../modules/local/braker3'
include { TSEBRA_TSEBRA as TSEBRA                               } from '../../../modules/local/tsebra/tsebra'
include { SAMTOOLS_MERGE                                        } from '../../../modules/local/samtools/merge'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    fasta: Path
    species: String
    orthodb_clade: String
    excluded_clades: List<String>
    excluded_species: List<String>
    training_proteins: Set<Path>
    short_read_sorted_bams_bais: Set<Record>?
    tsebra_gtfs: Set<Path>
    tsebra_hintsfiles: Set<Path>
}


workflow BRAKER {

    take:
    ch_input: Channel<Input>
    skip_orthodb_download: Boolean
    min_prot_db_seq_length: Integer

    main:
    
    // ----------------------------------------------------------
    // PREPARE PROTEIN TRAINING SET FOR BRAKER
    // ----------------------------------------------------------

    ch_proteins = TRAINING_PROTEIN_PREPARATION(
        ch_input,
        skip_orthodb_download,
        min_prot_db_seq_length
    )

    // remainder true means a left outer join
    ch_input = ch_input.join( ch_proteins, by: 'id', remainder: true )

    ch_has_mappings   = ch_input.filter{ rec -> rec.mappings != null }
    ch_leave_me_alone = ch_input.filter{ rec -> rec.mappings == null }
    
    ch_merge_me       = ch_has_mappings.filter{ rec -> rec.mappings.size() > 1 }
    ch_not_no_merge   = ch_has_mappings.filter{ rec -> rec.mappings.size() <= 1 }

    // ----------------------------------------------------------
    // MERGE MULTIPLE BAM FILES INTO A SINGLE BAM WHEN NECESSARY
    // ----------------------------------------------------------

    ch_samtools_merge_input = ch_merge_me.map{ rec ->
        def bams = rec.short_read_sorted_bams_bais.collect { r -> r.bam }
        def bais = rec.short_read_sorted_bams_bais.collect { r -> r.bai }
        record(id: rec.id, bams: bams, bais: bais)
    }

    ch_merged = SAMTOOLS_MERGE( ch_samtools_merge_input )

    ch_merge_me       = ch_merge_me.join( ch_merged, by: 'id' )
    ch_not_no_merge   = ch_not_no_merge.map { rec -> rec + record(bam: rec.short_read_sorted_bams_bais[0].bam) }

    ch_input = ch_leave_me_alone.mix( ch_not_no_merge ).mix( ch_merge_me )

    // ----------------------------------------------------------
    // RUN BRAKER3
    // ----------------------------------------------------------

    ch_braker_out = BRAKER3(
        ch_input.map { rec -> rec.subMap(['id', 'species', 'fasta', 'proteins_fasta', 'bam']) }
    )

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

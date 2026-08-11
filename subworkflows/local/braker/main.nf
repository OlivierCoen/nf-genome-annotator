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
    clade: String
    orthodb_clade: String
    excluded_clades: Iterable<String>
    excluded_species: Iterable<String>
    training_proteins: Iterable<Path>
    mappings: Iterable<Record>
    tsebra_gtfs: Iterable<Path>
    tsebra_hintsfiles: Iterable<Path>
}


workflow BRAKER {

    take:
    ch_input
    skip_orthodb_download
    min_prot_db_seq_length

    main:

    // ----------------------------------------------------------
    // PREPARE PROTEIN TRAINING SET FOR BRAKER
    // ----------------------------------------------------------

    TRAINING_PROTEIN_PREPARATION(
        ch_input.map { rec -> rec.subMap(['id', 'clade', 'orthodb_clade', 'excluded_clades', 'excluded_species', 'training_proteins']) },
        skip_orthodb_download,
        min_prot_db_seq_length
    )

    ch_input = ch_input.join( TRAINING_PROTEIN_PREPARATION.out.proteins, by: 'id' )

    // ----------------------------------------------------------
    // MERGE MULTIPLE BAM FILES INTO A SINGLE BAM WHEN NECESSARY
    // ----------------------------------------------------------

    ch_merge_me       = ch_input.filter{ rec -> rec.mappings.size() > 1 }
    ch_leave_me_alone = ch_input.filter{ rec -> rec.mappings.size() <= 1 }

    ch_samtools_merge_input = ch_merge_me.map{ rec ->
        def bams = rec.mappings.collect { r -> r.bam }
        def bais = rec.mappings.collect { r -> r.bai }
        record(id: rec.id, bams: bams, bais: bais)
    }

    SAMTOOLS_MERGE( ch_samtools_merge_input )

    ch_merge_me       = ch_merge_me.join( SAMTOOLS_MERGE.out, by: 'id' )
    ch_leave_me_alone = ch_leave_me_alone.map { rec -> rec + record(bam: null) }

    ch_input = ch_leave_me_alone.mix( ch_merge_me )

    // ----------------------------------------------------------
    // RUN BRAKER3
    // ----------------------------------------------------------

    BRAKER3(
        ch_input.map { rec -> rec.subMap(['id', 'species', 'fasta', 'training_proteins', 'bam']) }
    )

    ch_input = ch_input.join(BRAKER3.out, by: 'id')

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

    TSEBRA( ch_tsebra_input )

    ch_merged = ch_to_merge_with_tsebra.join( TSEBRA.out, by: 'id')

    // ----------------------------------------------------------
    // MIXING MERGED AND NOT MERGED
    // ----------------------------------------------------------

    ch_merged       = ch_merged.map      { rec -> rec + record(structural_annotation_gtf: rec.merged_gtf) }
    ch_not_to_merge = ch_not_to_merge.map{ rec -> rec + record(structural_annotation_gtf: rec.braker_gtf) }

    emit:
    annotated = ch_merged.mix( ch_not_to_merge )

}

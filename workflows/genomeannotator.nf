nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOME_PREPARATION                                            } from '../subworkflows/local/genome_preparation'
include { TAXONOMY_INFO                                                 } from '../subworkflows/local/taxonomy_info'
include { GENOME_MASKING                                                } from '../subworkflows/local/genome_masking'
include { PREPARE_RNASEQ_DATA                                           } from '../subworkflows/local/prepare_rnaseq_data'
include { STRUCTURAL_ANNOTATION                                         } from '../subworkflows/local/structural_annotation'
include { POST_PROCESS_ANNOTATION                                       } from '../subworkflows/local/post_process_annotation'
include { EXTRACT_SEQUENCES                                             } from '../subworkflows/local/extract_sequences'
include { FUNCTIONAL_ANNOTATION                                         } from '../subworkflows/local/functional_annotation'
include { QUALITY_CONTROLS                                              } from '../subworkflows/local/quality_controls'
include { REPORTING                                                     } from '../subworkflows/local/reporting'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Samplesheet {
    id: String
    fasta: Path
    species: String
    gff: Path?
    supplied_short_read_bams: Set<Path>
    supplied_short_reads: Set<Record>
    supplied_long_reads: Set<Record>
    supplied_short_read_sra_ids: List<String>
    supplied_long_read_sra_ids: List<String>
    training_proteins: Set<Path>
    orthodb_excluded_clades: List<String>
    orthodb_excluded_species: List<String>
    mmseqs_db: String
    tsebra_gtfs: Set<Path>
    tsebra_hintsfiles: Set<Path>
}


workflow GENOMEANNOTATOR {

    take:
    ch_main: Channel<Samplesheet>

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GENOME PREPARATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_prepared = GENOME_PREPARATION( ch_main )
    ch_main = ch_main.join( ch_prepared, by: 'id' )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FETCH NCBI TAXON ID, BUSCO DATASET AND ORTHODB CLADE
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_taxonomy = TAXONOMY_INFO(
        ch_main.map{ rec -> rec.species }.unique()
    )
    ch_main = ch_main.join( ch_taxonomy, by: 'species' )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GENOME MASKING
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !params.skip_masking ) {
        ch_masked = GENOME_MASKING (
            ch_main,
            params.genome_masker,
            params.dfam_db
        )
        ch_main = ch_main.join( ch_masked, by: 'id' )
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PREPARATION (QC / CLEANING / MAPPING) OF RNASEQ DATA (SHORT / LONG READS)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // only a subset of structural annotator can use RNAseq data 
    def structural_annotator_needing_rnaseq = ['braker3', 'tiberius']
    if ( params.structural_annotator in structural_annotator_needing_rnaseq || params.force_prepare_rnaseq ){

        ch_prepared_rnaseq_data = PREPARE_RNASEQ_DATA( 
            ch_main,
            params.skip_fastqc,
            params.skip_fastqc_raw,
            params.skip_fastqc_cleaned,
            params.skip_umi_extract,
            params.skip_short_read_cleaning,
            params.nb_short_read_sra_datasets,
            params.nb_long_read_sra_datasets,
            params.sra_max_size,
            params.sra_allow_single_end,
            params.sra_random_seed,
            params.read_sampling_size,
            params.short_read_mapper,
            params.ignore_existing_gff_for_mapping
        )
        ch_main = ch_main.join( ch_prepared_rnaseq_data, by: 'id', remainder: true )

    }
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STRUCTURAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !params.skip_structural_annotation ) {
    
        ch_structural_annotation = STRUCTURAL_ANNOTATION (
            ch_main,
            params.structural_annotator,
            params.mmseqs_db,
            params.skip_orthodb_download,
            params.skip_mmseqs_db_download,
            params.min_prot_db_seq_length
        )

        // NOTE: in case the structural annotation is performed,
        // samples for which annotation could not be performed are not kept for the following steps
        ch_main = ch_main.join( ch_structural_annotation, by: 'id' )

    } else {
        // IMPORTANT: when skipping the structural annotation, the provided reference gff becomes the structural annotation
        ch_main = ch_main.map { rec -> rec + record(structural_annotation: rec.reference_gff) }
    }

    // At this step, structural annotations are under the 'structural_annotation' key

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // POST-PROCESS STRUCTURAL ANNOTATION
    // COMPLEMENT WITH PREVIOUS ANNOTATION, CLEAN, MAKE ALTERNATIVE ANNOTATIONS, ...
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_post_processed_annotation = POST_PROCESS_ANNOTATION(
        ch_main.filter { rec -> rec.structural_annotation != null },
        params.complementation_use_ref_as_template,
        params.complement_annotation,
        params.skip_gff_cleaning,
        params.skip_alternative_annotations,
        params.gff_fix_feature_locations_duplicated,
        params.gff_fix_overlapping_genes,
        params.gff_filter_incomplete_gene_models
    )
    ch_main = ch_main.join( ch_post_processed_annotation, by: 'id', remainder: true )

    // At this step, cleaned structural annotations are under the 'gff' key

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE PROTEOME 
    // (ONLY IF FUNCTIONAL ANNOTATION OR QUALITY CONTROLS ARE TO BE RUN)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !(params.skip_functional_annotation && params.skip_qc) ) {
    
        ch_extracted_sequences = EXTRACT_SEQUENCES(
            ch_main.filter { rec -> rec.gff != null }
        )
        ch_main = ch_main.join( ch_extracted_sequences, by: 'id' )

    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FUNCTIONAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !params.skip_functional_annotation ) {

        // storing the final structural annotation
        ch_main = ch_main.map { rec -> rec + record(final_structural_annotation: rec.gff) }

        ch_functional_annotation = FUNCTIONAL_ANNOTATION (
            ch_main.filter { rec -> rec.proteome != null },
            params.functional_annotators,
            params.eggnog_mapper_mode,
            params.interproscan5_db,
            params.interproscan5_db_url
        )

        ch_main = ch_main.join( ch_functional_annotation, by: 'id' )

    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // VARIOUS QUALITY CONTROLS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !params.skip_qc ) {

        ch_qc = QUALITY_CONTROLS(
            ch_main,
            params.skip_busco,
            params.skip_omark,
            params.omamer_db_url,
            params.omamer_db
        )

    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MULTIQC & OTHER REPORTING
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_reporting = REPORTING(
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir
    )

    ch_main = ch_main.join( ch_reporting, by: 'id' )

    emit:
    results            = ch_main
    additional_results = channel.topic('additional_results')
    logs               = channel.topic('logs')

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

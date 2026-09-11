nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOME_PREPARATION                                            } from '../subworkflows/local/genome_preparation'
include { TAXONOMY_INFO                                                 } from '../subworkflows/local/taxonomy_info'
include { GENOME_MASKING                                                } from '../subworkflows/local/genome_masking'
include { FETCH_SRA_IDS                                                 } from '../subworkflows/local/fetch_sra_ids'
include { DOWNLOAD_READS                                                } from '../subworkflows/local/download_reads'
include { MAP_RNASEQ_READS                                              } from '../subworkflows/local/map_rnaseq_reads'
include { BAM_SORT_INDEX_STATS                                          } from '../subworkflows/local/bam_sort_index_stats'
include { STRUCTURAL_ANNOTATION                                         } from '../subworkflows/local/structural_annotation'
include { COMPLEMENT_ANNOTATION                                         } from '../subworkflows/local/complement_annotation'
include { CLEAN_ANNOTATION                                              } from '../subworkflows/local/clean_annotation'
include { ALTERNATIVE_ANNOTATIONS                                       } from '../subworkflows/local/alternative_annotation'
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
    supplied_rnaseq_bams: Iterable<Path>
    supplied_rnaseq_fastqs: Iterable<Record>
    rnaseq_experiment_ids: Iterable<String>
    training_proteins: Iterable<Path>
    orthodb_excluded_clades: Iterable<String>
    orthodb_excluded_species: Iterable<String>
    mmseqs_db: String
    tsebra_gtfs: Iterable<Path>
    tsebra_hintsfiles: Iterable<Path>
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

    if ( !params.skip_structural_annotation ) {

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
        // WHEN NEEDED, DOWNLOAD READS FROM PUBLIC DATABASES AND MAP THEM TO THE GENOME
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
        // only a subset of structural annotator can use RNAseq data 
        if ( params.structural_annotator in ['braker3', 'tiberius'] ){

            if ( params.fetch_sra_rnaseq ) {

                ch_sra_ids = FETCH_SRA_IDS( 
                    ch_main,
                    params.nb_short_read_sra_datasets,
                    params.nb_long_read_sra_datasets,
                    params.sra_allow_single_end,
                    params.sra_random_seed
                )
               // ch_main = ch_main.join( ch_sra_ids, by: 'id', remainder: true )
                
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // DOWNLOAD READS FROM SRA / ENA
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
            ch_downloaded_reads = DOWNLOAD_READS( ch_main ) 
            ch_main = ch_main.join( ch_downloaded_reads, by: 'id', remainder: true )
    
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // MAP RNASEQ READS TO GENOME
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
            ch_main = ch_main.map{ rec ->
                def downloaded_rnaseq_fastqs = rec.downloaded_rnaseq_fastqs ?: []
                rec + record(reads_to_map: rec.supplied_rnaseq_fastqs + downloaded_rnaseq_fastqs)
            }

            ch_reads_mapped = MAP_RNASEQ_READS(
                ch_main.filter{ rec -> rec.reads_to_map.size() > 0 }, // pass only samples for which there are reads
                params.skip_fastqc,
                params.skip_umi_extract,
                params.skip_trimming,
                params.rnaseq_mapper,
                params.ignore_existing_gff_for_mapping
            )
            ch_main = ch_main.join( ch_reads_mapped, by: 'id', remainder: true )

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // SORT ALL BAMS (SUPPLIED + NEWLY PRODUCED) AND GET MAPPING STATS
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
            ch_main = ch_main.map{ rec ->
                def new_rnaseq_bams = rec.new_rnaseq_bams ?: []
                rec + record(bams: rec.supplied_rnaseq_bams + new_rnaseq_bams)
            }

            ch_sorted_bam = BAM_SORT_INDEX_STATS(
                ch_main.filter { rec -> rec.bams.size() > 0 }
            )
            ch_main = ch_main.join( ch_sorted_bam, by: 'id', remainder: true )
            
        }
        
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STRUCTURAL ANNOTATION
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_structural_annotation = STRUCTURAL_ANNOTATION (
            ch_main,
            params.structural_annotator,
            params.mmseqs_db,
            params.skip_orthodb_download,
            params.skip_mmseqs_db_download,
            params.min_prot_db_seq_length
        )

        // NOTE: in case the structural annotation was performed
        // samples for which annotation could not be performed
        // are not kept for the following steps
        ch_main = ch_main.join( ch_structural_annotation, by: 'id')

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // COMPLEMENTATION OF ANNOTATION (WHEN NECESSARY)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        if ( params.complement_annotation ) {
            ch_complemented = COMPLEMENT_ANNOTATION( ch_main )
            ch_main = ch_main.join( ch_complemented, by: 'id' )
        }

    } else {
        // when skipping the structural annotation, the provided gff becomes the structural annotation
        ch_main = ch_main.map { rec -> rec + record(structural_annotation: rec.gff)}
    }

    // storing the provided gff (if any)
    // filtering to keep only records that have at least a structural annotation or a gff
    ch_main = ch_main
                .filter { rec -> rec.structural_annotation != null }
                .map { rec -> rec.gff ? rec + record(previous_annotation: rec.gff) : rec }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CLEANING OF GFF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !params.skip_gff_cleaning ) {
    
        ch_cleaned = CLEAN_ANNOTATION (
            ch_main,
            params.gff_fix_feature_locations_duplicated,
            params.gff_fix_overlapping_genes,
            params.gff_filter_incomplete_gene_models
        )
        ch_main = ch_main.join( ch_cleaned, by: 'id' )
    
        // NOTE: now the annotation is under the 'gff' key

    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE ALTERNATIVE ANNOTATIONS (LONGEST ISOFORMS ONLY, ...)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !params.skip_alternative_annotations ) {

        ch_alternative_annotations = ALTERNATIVE_ANNOTATIONS( ch_main )
        ch_main = ch_main.join( ch_alternative_annotations, by: 'id' )

    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE PROTEOME 
    // (ONLY IF FUNCTIONAL ANNOTATION OR QUALITY CONTROLS ARE TO BE RUN)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !(params.skip_functional_annotation && params.skip_qc) ) {
    
        ch_extracted_sequences = EXTRACT_SEQUENCES (
            ch_main
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
            ch_main,
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
    results = ch_main

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

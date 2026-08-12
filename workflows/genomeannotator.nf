nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOME_PREPARATION                                            } from '../subworkflows/local/genome_preparation'
include { TAXONOMY_INFO                                                 } from '../subworkflows/local/taxonomy_info'
include { GENOME_MASKING                                                } from '../subworkflows/local/genome_masking'
include { DOWNLOAD_READS                                                } from '../subworkflows/local/download_reads'
include { MAP_RNASEQ_READS                                              } from '../subworkflows/local/map_rnaseq_reads'
include { BAM_SORT_INDEX_STATS                                          } from '../subworkflows/local/bam_sort_index_stats'
include { STRUCTURAL_ANNOTATION                                         } from '../subworkflows/local/structural_annotation'
include { COMPLEMENT_ANNOTATION                                         } from '../subworkflows/local/complement_annotation'
include { CLEAN_ANNOTATIONS                                             } from '../subworkflows/local/clean_annotations'
include { ALTERNATIVE_ANNOTATIONS                                       } from '../subworkflows/local/alternative_annotation'
include { GET_PROTEOMES                                                 } from '../subworkflows/local/get_proteomes'
include { FUNCTIONAL_ANNOTATION                                         } from '../subworkflows/local/functional_annotation'
include { QUALITY_CONTROLS                                              } from '../subworkflows/local/qc'
include { REPORTING                                                     } from '../subworkflows/local/reporting'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Samplesheet {
    meta: Map
    genome: Path
}

workflow GENOMEANNOTATOR {

    take:
    ch_main: Channel<Samplesheet>

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GENOME PREPARATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_out = GENOME_PREPARATION( ch_main )
    ch_main = ch_main.join( ch_out.prepared, by: 'id' )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FETCH NCBI TAXON ID, BUSCO DATASET AND ORTHODB CLADE
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_out = TAXONOMY_INFO(
        ch_main.map{ rec -> rec.species }.unique()
    )
    ch_main = ch_main.join( ch_out.taxonomy, by: 'species' )

    if ( !params.skip_structural_annotation ) {

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // GENOME MASKING
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if ( !params.skip_masking ) {
            ch_out = GENOME_MASKING (
                ch_main,
                params.genome_masker
            )
            ch_main = ch_main.join( ch_out.masked, by: 'id' )
        }

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // DOWNLOAD READS FROM SRA / ENA IF NEEDED
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_out = DOWNLOAD_READS( ch_main ) 
        ch_main = ch_main.join( ch_out.reads, by: 'id' )

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // MAP RNASEQ READS TO GENOME
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_main = ch_main.map{ rec ->
            def downloaded_rnaseq_fastqs = rec.downloaded_rnaseq_fastqs ?: []
            rec + record(reads_to_map: rec.supplied_rnaseq_fastqs + downloaded_rnaseq_fastqs)
        }

        ch_out = MAP_RNASEQ_READS(
            ch_main.filter{ rec -> rec.reads_to_map.size() > 0 }, // pass only samples for which there are reads
            params.skip_fastqc,
            params.skip_umi_extract,
            params.skip_trimming,
            params.rnaseq_mapper,
            params.ignore_existing_gff_for_mapping
        )
        ch_main = ch_main.join( ch_out.mapped, by: 'id' )

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // SORT ALL BAMS (SUPPLIED + NEWLY PRODUCED) AND GET MAPPING STATS
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_main = ch_main.map{ rec ->
            def new_rnaseq_bams = rec.new_rnaseq_bams ?: []
            rec + record(bams: rec.supplied_rnaseq_bams + new_rnaseq_bams)
        }

        ch_out = BAM_SORT_INDEX_STATS(
            ch_main.map { rec -> rec.bams.size() > 0 }
        )
        ch_main = ch_main.join( ch_out.sorted_indexed, by: 'id' )

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STRUCTURAL ANNOTATION
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_out = STRUCTURAL_ANNOTATION (
            ch_main,
            params.structural_annotator,
            params.mmseqs_db,
            params.skip_orthodb_download,
            params.skip_mmseqs_db_download,
            params.min_prot_db_seq_length
        )
        ch_main = ch_main.join( ch_out.annotated, by: 'id')

    } else {
        // when skipping the structural annotation, the provided gff becomes the structural annotation
        ch_main = ch_main.map { rec -> rec + record(structural_annotation: rec.gff)}
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // COMPLEMENTATION OF ANNOTATION (WHEN NECESSARY)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( params.complement_annotation ) {
        ch_out = COMPLEMENT_ANNOTATION( ch_main )
        ch_main = ch_main.join( ch_out.complemented, by: 'id' )
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CLEANING OF GTF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*
    CLEAN_ANNOTATIONS (
        ch_structural_annotations,
        ch_genome,
        params.gff_fix_feature_locations_duplicated,
        params.skip_gff_fix_overlapping_genes,
        params.skip_gff_filter_incomplete_gene_models
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE ALTERNATIVE ANNOTATIONS (LONGEST ISOFORMS ONLY, ...)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ALTERNATIVE_ANNOTATIONS( ch_gff )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ORGANISE ANNOTATIONS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_structural_annotation = CLEAN_ANNOTATIONS.out.gff
                            .map {
                                meta, file -> [ meta + [final_annotation: true], file ]
                            }

    // for now, the final annotation is set to the structural annotation
    // is will be set to the functional annotation if it is not skipped
    ch_final_annotation = ch_structural_annotation

    ch_intermediate_annotations = ch_structural_annotations
                                    .mix( CLEAN_ANNOTATIONS.out.intermediate_gffs )
                                    .map {
                                        meta, file -> [ meta + [final_annotation: false], file ]
                                    }

    ch_alternative_annotations = ALTERNATIVE_ANNOTATIONS.out.longest_isoforms_gff
                                    .map {
                                        meta, file -> [ meta + [final_annotation: false], file ]
                                    }

    ch_all_annotations = ch_structural_annotation
                            .mix( ch_intermediate_annotations )
                            .mix( ch_alternative_annotations )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE PROTEOME
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GET_PROTEOMES (
        ch_all_annotations,
        ch_genome,
        params.codon_usage_id
    )

    ch_proteomes = GET_PROTEOMES.out.proteomes
    ch_main_proteome = ch_proteomes
                        .filter{ meta, file -> meta.final_annotation == true }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FUNCTIONAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !params.skip_functional_annotation ) {

        FUNCTIONAL_ANNOTATION (
            ch_main_proteome,
            ch_structural_annotation,
            params.functional_annotators,
            params.interproscan_db,
            params.interproscan_db_url
        )

        ch_functional_annotation = FUNCTIONAL_ANNOTATION.out.gff
        ch_final_annotation      = ch_functional_annotation

        ch_versions = ch_versions
                        .mix( FUNCTIONAL_ANNOTATION.out.versions )
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // VARIOUS QUALITY CONTROLS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    QUALITY_CONTROLS (
        ch_genome,
        ch_busco_lineage,
        ch_all_annotations,
        ch_main_proteome,
        ch_proteomes,
        ch_structural_annotation,
        ch_functional_annotation,
        params.skip_omark,
        params.omamer_db_url,
        params.omamer_db
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MULTIQC & OTHER REPORTING
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    REPORTING(
        ch_versions,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir
    )

*/

}
    emit:
    results = ch_main

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

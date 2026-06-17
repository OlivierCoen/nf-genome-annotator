/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AGAT_SPCOMPLEMENTANNOTATIONS as COMPLEMENT_ANNOTATIONS        } from '../modules/local/agat/spcomplementannotations'


include { GENOME_PREPARATION                                            } from '../subworkflows/local/genome_preparation'
include { GENOME_MASKING                                                } from '../subworkflows/local/genome_masking'
include { DOWNLOAD_READS                                                } from '../subworkflows/local/download_reads'
include { MAP_TO_GENOME_SORT_INDEX                                      } from '../subworkflows/local/map_to_genome_sort_index'
include { STRUCTURAL_ANNOTATION                                         } from '../subworkflows/local/structural_annotation'
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

workflow GENOME_ANNOTATOR {

    take:
    ch_samplesheet

    main:

    ch_input = ch_samplesheet.multiMap{
                    meta, genome, gff, rnaseq_bam, rnaseq_fastq, rnaseq_ids, proteins, braker_gtf, braker_hintsfile ->
                        genome: [ meta, genome ]
                        gff: gff ? [ meta, gff ] : [[:], []]
                        rnaseq_bam: rnaseq_bam ? [ meta, rnaseq_bam ] : [[:], []]
                        rnaseq_fastq: rnaseq_fastq ? [ meta, rnaseq_fastq ] : [[:], []]
                        rnaseq_id: rnaseq_ids ? [ meta, rnaseq_ids ] : [[:], []]
                        protein: proteins ? [ meta, proteins ] : [[:], []]
                        braker_gtf: braker_gtf ? [ meta, braker_gtf ] : [[:], []]
                        braker_hintsfile: braker_hintsfile ? [ meta, braker_hintsfile ] : [[:], []]
                }

    ch_genome       = ch_input.genome

    ch_gff          = ch_input.gff
                        .filter { meta, file -> file != []}

    ch_proteins     = ch_input.protein
                        .transpose()
                        .filter { meta, fasta -> fasta != [] }
                        .groupTuple()

    ch_provided_rnaseq_fastq = ch_input.rnaseq_fastq
                                .transpose()
                                .filter { meta, reads -> reads != []}
                                .map { meta, reads ->
                                    fastq_1 = reads[0]
                                    fastq_2 = reads[1]
                                    if ( fastq_2 ) {
                                        [ meta + [ single_end: false ], [ fastq_1, fastq_2 ] ]
                                    } else {
                                        [ meta + [ single_end: true ], fastq_1 ]
                                    }
                                }

    ch_rnaseq_bam   = ch_input.rnaseq_bam
                        .transpose()
                        .filter { meta, file -> file != []}

    ch_braker_gtf   = ch_input.braker_gtf
                        .filter { meta, file -> file != []}

    ch_braker_hintsfile    = ch_input.braker_hintsfile
                        .filter { meta, file -> file != []}

    ch_rnaseq_id   = ch_input.rnaseq_id
                        .transpose()
                        .filter { meta, id -> id != []}


    ch_versions = channel.empty()

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GENOME PREPARATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GENOME_PREPARATION ( ch_genome )
    ch_genome = GENOME_PREPARATION.out.prepared_genome


    if ( !params.skip_structural_annotation ) {

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // GENOME MASKING
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if ( !params.skip_masking ) {
            GENOME_MASKING ( ch_genome )
            ch_genome = GENOME_MASKING.out.masked_genome
        }

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // DOWNLOAD READS FROM SRA / ENA IF NEEDED
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        DOWNLOAD_READS( ch_rnaseq_id )

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // MAP RNASEQ READS TO GENOME
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_rnaseq_fastq = ch_provided_rnaseq_fastq
                            .mix( DOWNLOAD_READS.out.reads )

        MAP_TO_GENOME_SORT_INDEX(
            ch_genome,
            ch_rnaseq_fastq,
            ch_rnaseq_bam,
            ch_gff,
            params.skip_fastqc,
            params.skip_umi_extract,
            params.skip_trimming,
            params.rnaseq_mapper,
            params.ignore_existing_gtf_for_mapping
        )

        ch_grouped_bam_bai = MAP_TO_GENOME_SORT_INDEX.out.bam_bai
                                .groupTuple()

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STRUCTURAL ANNOTATION
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        STRUCTURAL_ANNOTATION (
            ch_genome,
            ch_proteins,
            ch_grouped_bam_bai,
            ch_braker_gtf,
            ch_braker_hintsfile,
            params.structural_annotator,
            params.species,
            params.busco_lineage,
            params.clade,
            params.excluded_clades,
            params.excluded_species,
            params.mmseqs_db,
            params.skip_orthodb_download,
            params.skip_mmseqs_db_download,
            params.min_prot_db_seq_length
        )

        ch_structural_annotations = STRUCTURAL_ANNOTATION.out.annotations

        ch_versions = ch_versions
                        .mix( STRUCTURAL_ANNOTATION.out.versions )

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // COMPLEMENTATION OF ANNOTATION (WHEN NECESSARY)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if ( !params.complement_annotation ) {

            ch_branched_annotations = ch_structural_annotations
                                        .join( ch_gff, remainder: true )
                                        .branch{
                                            meta, annotation, gff ->
                                                to_complement: gff != null
                                                    [ meta, gff, annotation ]
                                                leave_me_alone: gff == null
                                                    [ meta, annotation ]
                                        }

            COMPLEMENT_ANNOTATIONS ( ch_branched_annotations.to_complement, [] )

            ch_annotation = ch_branched_annotations.leave_me_alone
                                .mix( COMPLEMENT_ANNOTATIONS.out.gff )

        }

    } else {
        ch_structural_annotations = ch_gff
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CLEANING OF GTF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    ch_gff = CLEAN_ANNOTATIONS.out.gff
                .map {
                    meta, file -> [ meta + [final_annotation: true], file ]
                }

    ch_intermediate_annotations = ch_structural_annotations
                                    .mix( CLEAN_ANNOTATIONS.out.intermediate_gffs )
                                    .map {
                                        meta, file -> [ meta + [final_annotation: false], file ]
                                    }

    ch_alternative_annotations = ALTERNATIVE_ANNOTATIONS.out.longest_isoforms_gff
                                    .map {
                                        meta, file -> [ meta + [final_annotation: false], file ]
                                    }

    ch_all_annotations = ch_gff
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
            ch_gff,
            params.functional_annotators,
            params.interproscan_db,
            params.interproscan_db_url
        )

        ch_functional_annotation = FUNCTIONAL_ANNOTATION.out.gff

        ch_versions = ch_versions
                        .mix( FUNCTIONAL_ANNOTATION.out.versions )
    } else {
        ch_functional_annotation = channel.empty()
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // VARIOUS QUALITY CONTROLS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    QUALITY_CONTROLS (
        ch_genome,
        ch_all_annotations,
        ch_main_proteome,
        ch_proteomes,
        ch_gff,
        ch_functional_annotation,
        params.skip_omark,
        params.omamer_db_url
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


    emit:
    multiqc_report = REPORTING.out.report.toList()
    versions       = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AGAT_SPCOMPLEMENTANNOTATIONS as COMPLEMENT_ANNOTATIONS        } from '../modules/local/agat/spcomplementannotations'


include { GENOME_PREPARATION                                            } from '../subworkflows/local/genome_preparation'
include { GENOME_MASKING                                                } from '../subworkflows/local/genome_masking'
include { MAP_TO_GENOME_SORT_INDEX                                      } from '../subworkflows/local/map_to_genome_sort_index'
include { STRUCTURAL_ANNOTATION                                         } from '../subworkflows/local/structural_annotation'
include { CLEAN_ANNOTATIONS                                             } from '../subworkflows/local/clean_annotations'
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
                    meta, genome, rnaseq_bam, rnaseq_fastq, rnaseq_sra, proteins, gtf, hintsfile, ref_gff ->
                        genome: [ meta, genome ]
                        rnaseq_bam: rnaseq_bam ? [ meta, rnaseq_bam ] : [[:], []]
                        rnaseq_fastq: rnaseq_fastq ? [ meta, rnaseq_fastq ] : [[:], []]
                        rnaseq_sra: rnaseq_sra ? [ meta, rnaseq_sra ] : [[:], []]
                        proteins: proteins ? [ meta, proteins ] : [[:], []]
                        gtf: gtf ? [ meta, gtf ] : [[:], []]
                        hintsfile: hintsfile ? [ meta, hintsfile ] : [[:], []]
                        ref_gff: ref_gff ? [ meta, ref_gff ] : [[:], []]
                }

    ch_genome       = ch_input.genome

    ch_proteins     = ch_input.proteins
                        .transpose()
                        .filter { meta, fasta -> fasta != [] }
                        .groupTuple()

    ch_rnaseq_fastq = ch_input.rnaseq_fastq
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
                        .filter { meta, bam -> bam != []}

    ch_gtf          = ch_input.gtf
                        .filter { meta, gtf -> gtf != []}

    ch_hintsfile    = ch_input.hintsfile
                        .filter { meta, file -> file != []}

    //ch_ref_gff      = ch_input.ref_gff.filter { meta, gff -> gff != []}.ifEmpty([])

    ch_rnaseq_sra   = ch_input.rnaseq_sra
                        .filter { meta, sra -> sra != []}


    ch_versions = channel.empty()

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GENOME PREPARATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GENOME_PREPARATION ( ch_genome )
    ch_genome = GENOME_PREPARATION.out.prepared_genome

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GENOME MASKING
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !params.skip_masking ) {
        GENOME_MASKING ( ch_genome )
        ch_genome = GENOME_MASKING.out.masked_genome
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP RNASEQ READS TO GENOME
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    MAP_TO_GENOME_SORT_INDEX(
        ch_genome,
        ch_rnaseq_fastq,
        ch_rnaseq_bam,
        ch_gtf,
        params.skip_fastqc,
        params.skip_umi_extract,
        params.skip_trimming,
        params.star_ignore_existing_gtf
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
        ch_gtf,
        ch_hintsfile,
        params.structural_annotator,
        params.species,
        params.clade,
        params.excluded_clades,
        params.excluded_species
    )
    ch_structural_annotations = STRUCTURAL_ANNOTATION.out.annotations

    /*
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // COMPLEMENTATION OF ANNOTATION (WHEN NECESSARY)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_branched_annotations = ch_structural_annotations
                                .join( ch_ref_gff, remainder: true )

                                .view{ v -> "after join $v"}
                                .branch{
                                    meta, annotation, ref_gff ->
                                        to_complement: ref_gff != null
                                            [ meta, ref_gff, annotation ]
                                        leave_me_alone: ref_gff == null
                                            [ meta, annotation ]
                                }

    COMPLEMENT_ANNOTATIONS ( ch_branched_annotations.to_complement, [] )

    ch_annotation = ch_branched_annotations.leave_me_alone
                        .mix( COMPLEMENT_ANNOTATIONS.out.gff )
    */
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CLEANING OF GTF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    CLEAN_ANNOTATIONS (
        ch_structural_annotations,
        ch_genome,
        params.skip_gff_keep_longest_isoform,
        params.skip_gff_fix_overlapping_genes,
        params.skip_gff_filter_incomplete_gene_models,
        params.skip_gff_fix_cds_phases
    )

    ch_gff = CLEAN_ANNOTATIONS.out.gff
                .map {
                    meta, file -> [ meta + [final_annotation: true], file ]
                }

    ch_intermediate_annotations = ch_structural_annotations
                                    .mix( CLEAN_ANNOTATIONS.out.intermediate_gffs )
                                    .map {
                                        meta, file -> [ meta + [final_annotation: false], file ]
                                    }

    ch_all_annotations = ch_gff.mix( ch_intermediate_annotations )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE PROTEOME
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GET_PROTEOMES (
        ch_all_annotations,
        ch_genome,
        params.codon_usage_id
    )

    ch_proteomes = GET_PROTEOMES.out.proteomes
    ch_main_proteome = ch_proteomes.filter{ meta, file -> meta.final_annotation == true }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FUNCTIONAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !params.skip_functional_annotation ) {
        FUNCTIONAL_ANNOTATION (
            ch_main_proteome,
            ch_gff,
            params.functional_annotator,
            params.interproscan_db,
            params.interproscan_db_url
        )
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // QUALITY CONTROLS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    QUALITY_CONTROLS (
        ch_genome,
        ch_all_annotations,
        ch_proteomes
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MULTIQC
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_versions = ch_versions
                    .mix( STRUCTURAL_ANNOTATION.out.versions )
                    .mix( FUNCTIONAL_ANNOTATION.out.versions )

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

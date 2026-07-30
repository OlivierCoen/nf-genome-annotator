nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AGAT_SPCOMPLEMENTANNOTATIONS as COMPLEMENT_ANNOTATIONS        } from '../modules/local/agat/spcomplementannotations'


include { GENOME_PREPARATION                                            } from '../subworkflows/local/genome_preparation'
include { TAXONOMY_INFO                                                 } from '../subworkflows/local/taxonomy_info'
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

record Samplesheet {
    meta: Map
    genome: Path
}

workflow GENOMEANNOTATOR {

    take:
    ch_samplesheet: Channel<Samplesheet>

    main:

    ch_main = ch_samplesheet
                .map{ meta, genome ->
                    record(
                        id: meta.id,
                        fasta: genome,
                        species: meta.species,
                        gff: meta.gff ?: [],
                        rnaseq_bam: meta.rnaseq_bam ?: [],
                        rnaseq_fastq: meta.rnaseq_fastq ?: [],
                        rnaseq_public_id: meta.rnaseq_public_id ?: [],
                        proteins: meta.proteins ?: [],
                        braker_gtf: meta.braker_gtf ?: [],
                        hintsfile: meta.hintsfile ?: []
                    )
                }

                /*
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
                                    def fastq_1 = reads[0]
                                    def fastq_2 = reads[1]
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
*/

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GENOME PREPARATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GENOME_PREPARATION ( 
        ch_main.map{ rec -> record(id: rec.id, fasta: rec.fasta) }
    )
    ch_main = ch_main.join(GENOME_PREPARATION.out.prepared, by: 'id')
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FETCH NCBI TAXON ID, BUSCO DATASET AND ORTHODB CLADE
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    TAXONOMY_INFO( 
        ch_main.map{ rec -> rec.species }.unique()
    )
    ch_main = ch_main.join(TAXONOMY_INFO.out.taxonomy, by: 'species')

    if ( !params.skip_structural_annotation ) {

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // GENOME MASKING
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if ( !params.skip_masking ) {
            GENOME_MASKING ( 
                ch_main.map{ rec -> record(id: rec.id, fasta: rec.fasta) }
            )
            ch_main = ch_main.join(GENOME_MASKING.out.masked, by: 'id')
        }
/*
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // DOWNLOAD READS FROM SRA / ENA IF NEEDED
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        DOWNLOAD_READS( ch_input )
        ch_downloaded_reads = DOWNLOAD_READS.out.reads

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // MAP RNASEQ READS TO GENOME
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_rnaseq_fastq = ch_provided_rnaseq_fastq
                            .mix( ch_downloaded_reads )

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

        ch_bam_bai = MAP_TO_GENOME_SORT_INDEX.out.bam_bai
                               
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // STRUCTURAL ANNOTATION
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        STRUCTURAL_ANNOTATION (
            ch_genome,
            ch_proteins,
            ch_bam_bai.groupTuple(),
            ch_braker_gtf,
            ch_braker_hintsfile,
            params.structural_annotator,
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
*/
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

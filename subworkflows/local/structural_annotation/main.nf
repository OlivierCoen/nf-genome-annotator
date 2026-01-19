include { BUSCO_DOWNLOADPROTEINS as DOWNLOAD_ORTHODB_PROTEINS   } from '../../../modules/local/busco/download'
include { BRAKER3                                               } from '../../../modules/local/braker3'
include { TSEBRA_TSEBRA          as TSEBRA                      } from '../../../modules/local/tsebra/tsebra'

include { MMSEQS_DATABASES                                      } from '../../../modules/nf-core/mmseqs/databases'
include { METAEUK_EASYPREDICT                                   } from '../../../modules/local/metaeuk/easypredict'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow STRUCTURAL_ANNOTATION {

    take:
    ch_genome

    main:

    ch_versions = channel.empty()

    // ----------------------------------------------------------
    // EXTRACT GENOMES ALREADY ASSOCIATED WITH A CUSTOM PROTEOME FILE
    // ----------------------------------------------------------

    ch_separated_on_proteins = ch_genome
                                .branch {
                                    meta, genome ->
                                        no_proteins: meta.proteins == []
                                        with_proteins: meta.proteins != []
                                            [ meta, genome, meta.proteins ]
                                }

    // ----------------------------------------------------------
    // DOWNLOAD PROTEIN DATABASE IF NECESSARY
    // ----------------------------------------------------------

    if ( params.mmseqs_db == null ) {

        // trick to execute the process DOWNLOAD_ORTHODB_PROTEINS
        // only if there are elements in the channel ch_separated_on_proteins.no_proteins
        def orthodb_lineage = params.orthodb_lineage ?: params.busco_lineage
        ch_separated_on_proteins.no_proteins
            .take (1)
            .map { orthodb_lineage }
            | DOWNLOAD_ORTHODB_PROTEINS

        ch_protein_db = DOWNLOAD_ORTHODB_PROTEINS.out.proteins

    } else { // params.orthodb_lineage == null && params.mmseqs_db != null

        ch_separated_on_proteins.no_proteins
            .take (1)
            .map { params.mmseqs_db }
            | MMSEQS_DATABASES

        ch_protein_db = MMSEQS_DATABASES.out.database

        ch_versions = ch_versions.mix( MMSEQS_DATABASES.out.versions )

    }

    // ----------------------------------------------------------
    // ADD DOWNLOAD PROTEINS TO THE ITEMS WITHOUT PROTEINS
    // ----------------------------------------------------------

    ch_no_proteins_with_protein_db = ch_separated_on_proteins.no_proteins
                                        .combine( ch_protein_db )
                                        .map {
                                            meta, genome, protein_db ->
                                                [ meta, genome, protein_db ]
                                        }

    ch_to_annotate = ch_separated_on_proteins.with_proteins
                        .mix ( ch_no_proteins_with_protein_db )

    if ( params.structural_annotator == "braker3" ) {

        // ----------------------------------------------------------
        // RUN BRAKER3 & TSEBRA
        // ----------------------------------------------------------

        // extracts bam file from the meta map
        BRAKER3(
            ch_to_annotate.map { meta, genome, proteins -> [ meta, genome, meta.bam, proteins ] }
        )

        // MERGE ANNOTATIONS WHEN NECESSARY

        // separate inputs that need to be merged from the rest
        ch_braker_gtfs = BRAKER3.out.gtf
                            .branch {
                                meta, gtf ->
                                    to_merge: meta.gtf != [] && meta.hintsfile != []
                                    not_to_merge: meta.gtf == [] || meta.hintsfile == []
                            }

        // add the newly computed hintsfile from BRAKER to these inputs
        ch_tsebra_input = ch_braker_gtfs.to_merge
                            .join ( BRAKER3.out.hintsfile )
                            .map {
                                meta, gtf, hintsfile ->
                                    [ meta, [ gtf, meta.gtf ], [ hintsfile, meta.hintsfile ] ]
                            }

        // TSEBRA
        TSEBRA(
            ch_tsebra_input,
            [],
            []
        )

        ch_annotations = ch_braker_gtfs.not_to_merge
                                    .mix( TSEBRA.out.merged_gtf )

    } else {

        // ----------------------------------------------------------
        // RUN METAEUK
        // ----------------------------------------------------------

        METAEUK_EASYPREDICT( ch_to_annotate )
        ch_annotations = METAEUK_EASYPREDICT.out.gff

    }

    emit:
    annotations             = ch_annotations
    versions                = ch_versions
}

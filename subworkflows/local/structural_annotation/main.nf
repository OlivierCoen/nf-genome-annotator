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

    ch_versions = Channel.empty()

    // ----------------------------------------------------------
    // EXTRACT GENOMES ALREADY ASSOCIATED WITH A CUSTOM PROTEOME FILE
    // ----------------------------------------------------------

    ch_genome
        .branch {
            meta, genome ->
                no_proteins: meta.proteins == []
                with_proteins: meta.proteins != []
                    [ meta, genome, meta.proteins ]
        }
        .set { ch_separated_on_proteins }

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

        DOWNLOAD_ORTHODB_PROTEINS.out.proteins.set { ch_protein_db }

    } else { // params.orthodb_lineage == null && params.mmseqs_db != null

        ch_separated_on_proteins.no_proteins
            .take (1)
            .map { params.mmseqs_db }
            | MMSEQS_DATABASES

        MMSEQS_DATABASES.out.database.set { ch_protein_db }

        ch_versions = ch_versions.mix( MMSEQS_DATABASES.out.versions )

    }

    // ----------------------------------------------------------
    // ADD DOWNLOAD PROTEINS TO THE ITEMS WITHOUT PROTEINS
    // ----------------------------------------------------------

    ch_separated_on_proteins.no_proteins
        .combine( ch_protein_db )
        .map {
            meta, genome, protein_db ->
                [ meta, genome, protein_db ]
        }
        .set { ch_no_proteins_with_protein_db }

    ch_separated_on_proteins.with_proteins
        .mix ( ch_no_proteins_with_protein_db )
        .set { ch_to_annotate }

    // ----------------------------------------------------------
    // CHOOSE STRUCTURAL ANNOTATORS
    // ----------------------------------------------------------

    if ( params.structural_annotator == "auto" ) { // auto mode for structural annotator

        ch_to_annotate
            .branch {
                meta, genome, protein_db ->
                    def threshold = params.genome_max_size_auto_metaeuk.toLong()
                    def genome_size = meta.genome_size.toLong()
                    metaeuk: genome_size < threshold
                    braker: genome_size >= threshold
            }
            .set { ch_to_annotate }

    } else { // put all genomes in the same channel (depending on params.structural_annotator)

        ch_to_annotate
            .branch {
                meta, genome, protein_db ->
                    metaeuk: params.structural_annotator == "metaeuk"
                    braker: params.structural_annotator == "braker3"
            }
            .set { ch_to_annotate }

    }

    // ----------------------------------------------------------
    // RUN BRAKER3 & TSEBRA
    // ----------------------------------------------------------

    // BRAKER3
    // extracts bam file from the meta map
    BRAKER3(
        ch_to_annotate.braker.map { meta, genome, proteins -> [ meta, genome, meta.bam, proteins ] }
    )

    // MERGE ANNOTATIONS WHEN NECESSARY

    // separate inputs that need to be merged from the rest
    BRAKER3.out.gtf
        .branch {
            meta, gtf ->
                to_merge: meta.gtf != [] && meta.hintsfile != []
                not_to_merge: meta.gtf == [] || meta.hintsfile == []
        }
        .set { ch_braker_gtfs }

    // add the newly computed hintsfile from BRAKER to these inputs
    ch_braker_gtfs.to_merge
        .join ( BRAKER3.out.hintsfile )
        .map {
            meta, gtf, hintsfile ->
                [ meta, [ gtf, meta.gtf ], [ hintsfile, meta.hintsfile ] ]
        }
        .set { ch_tsebra_input }

    // TSEBRA
    TSEBRA( ch_tsebra_input, [], [] )

    ch_braker_gtfs.not_to_merge
        .mix( TSEBRA.out.merged_gtf )
        .set { ch_braker_annotations }

    // ----------------------------------------------------------
    // RUN METAEUK
    // ----------------------------------------------------------

    METAEUK_EASYPREDICT( ch_to_annotate.metaeuk )
    METAEUK_EASYPREDICT.out.gff.set { ch_metaeuk_annotations }


    // ----------------------------------------------------------
    // MERGE ALL ANNOTATIONS
    // ----------------------------------------------------------

    ch_braker_annotations
        .mix( ch_metaeuk_annotations )
        .set { ch_annotations }

    emit:
    annotations             = ch_annotations
    versions                = ch_versions
}

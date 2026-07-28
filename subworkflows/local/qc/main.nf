include { BUSCO_DOWNLOAD                                              } from '../../../modules/local/busco/download'
include { BUSCO_BUSCO as BUSCO_GENOME                                 } from '../../../modules/local/busco/busco'
include { BUSCO_BUSCO as BUSCO_PROTEOME                               } from '../../../modules/local/busco/busco'
include { AGAT_SPSTATISTICS as AGAT_GTF_STATISTICS                    } from '../../../modules/local/agat/spstatistics'
include { AGAT_SPFUNCTIONALSTATISTICS as AGAT_FUNCTIONAL_STATISTICS   } from '../../../modules/local/agat/spfunctionalstatistics'

include { OMARK                                                       } from '../omark'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow QUALITY_CONTROLS {

    take:
    ch_genome
    ch_busco_lineage
    ch_all_annotations
    ch_main_proteome
    ch_all_proteomes
    ch_gff
    ch_functional_annotation
    skip_omark
    omamer_db_url
    omamer_db

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DOWNLOAD NECESSARY BUSCO DATASETS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    BUSCO_DOWNLOAD( 
        ch_busco_lineage.map{ meta, lineage -> lineage }.unique()
    )

    ch_busco_download = ch_busco_lineage
                            .combine( BUSCO_DOWNLOAD.out.download_dir )
                            .filter { meta, lineage1, lineage2, busco_downloads -> lineage1 == lineage2 }
                            .map { meta, lineage1, lineage2, busco_downloads -> [ meta, busco_downloads ] }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // BUSCO
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    BUSCO_GENOME (
        ch_genome.join( ch_busco_download ),
        'genome'
    )

    BUSCO_PROTEOME (
        ch_all_proteomes.join( ch_busco_download ),
        'proteins'
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ASSESSMENT OF ANNOTATION QUALITY WITH OMARK
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !skip_omark ) {

        OMARK(
            ch_main_proteome,
            ch_gff,
            omamer_db_url,
            omamer_db
        )

    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // METRICS OF STRUCTURAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    AGAT_GTF_STATISTICS ( ch_all_annotations )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // METRICS OF FUNCTIONAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    AGAT_FUNCTIONAL_STATISTICS( ch_functional_annotation )

}

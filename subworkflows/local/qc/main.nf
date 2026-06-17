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
    ch_all_annotations
    ch_main_proteome
    ch_all_proteomes
    ch_gff
    ch_functional_annotation
    skip_omark
    omamer_db_url

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // BUSCO
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def busco_lineages_path = []
    def busco_config_file = []
    def busco_clean_intermediates = true

    BUSCO_GENOME (
        ch_genome,
        'genome',
        params.busco_lineage,
        busco_lineages_path,
        busco_config_file,
        busco_clean_intermediates
    )

    BUSCO_PROTEOME (
        ch_all_proteomes,
        'proteins',
        params.busco_lineage,
        busco_lineages_path,
        busco_config_file,
        busco_clean_intermediates
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ASSESSMENT OF ANNOTATION QUALITY WITH OMARK
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !skip_omark ) {

        OMARK(
            ch_main_proteome,
            ch_gff,
            omamer_db_url
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

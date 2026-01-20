include { BUSCO_BUSCO as BUSCO_GENOME                       } from '../../../modules/local/busco/busco'
include { BUSCO_BUSCO as BUSCO_PROTEOME                     } from '../../../modules/local/busco/busco'
include { AGAT_SPSTATISTICS as AGAT_GTF_STATISTICS          } from '../../../modules/local/agat/spstatistics'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow QUALITY_CONTROLS {

    take:
    ch_genome
    ch_gtf
    ch_proteome

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
        ch_proteome,
        'proteins',
        params.busco_lineage,
        busco_lineages_path,
        busco_config_file,
        busco_clean_intermediates
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // QUALITY CONTROLS OF GTF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    AGAT_GTF_STATISTICS ( ch_gtf )

}

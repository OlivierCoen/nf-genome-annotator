include { AGAT_SPSTATISTICS as GTF_STATISTICS          } from '../../../modules/local/agat/spstatistics'
include { BUSCO_BUSCO as BUSCO_GENOME                  } from '../../../modules/local/busco/busco'
include { BUSCO_BUSCO as BUSCO_PROTEOME                } from '../../../modules/local/busco/busco'



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

    BUSCO_GENOME.out.short_summaries_txt
        .mix ( BUSCO_PROTEOME.out.short_summaries_txt )
        .set { ch_busco_short_summaries }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // QUALITY CONTROLS OF GTF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GTF_STATISTICS ( ch_gtf )



    emit:
    gtf_stats                           = GTF_STATISTICS.out.stats_yaml
    busco_short_summaries               = ch_busco_short_summaries

}

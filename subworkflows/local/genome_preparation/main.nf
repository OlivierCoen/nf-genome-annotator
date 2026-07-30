include { CHECK_GENOME              } from '../../../modules/local/check/genome'
include { SEQKIT_STATS              } from '../../../modules/local/seqkit/stats'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow GENOME_PREPARATION {

    take:
    ch_main

    main:

    // ----------------------------------------------------------
    // CHECK HEADERS OR WHOLE PROTEIN DB AND RAISE ERROR IF UNHANDLED CHARACTERS
    // ----------------------------------------------------------

    CHECK_GENOME( ch_main )
    ch_main = ch_main.join(CHECK_GENOME.out.checked, by: 'id')
    
    // ----------------------------------------------------------
    // COMPUTE STATS ABOUT GENOME
    // ----------------------------------------------------------

    SEQKIT_STATS ( ch_main )
    ch_stats = SEQKIT_STATS.out.stats

    ch_main = ch_main
                .join(ch_stats, by: 'id')
                .map {
                    rec ->
                        def csv = rec.genome_stats.splitCsv( header: true, sep: '\t', limit: 1 ).collect()
                        rec + record(genome_size: csv[0].sum_len.toInteger())
                }

    emit:
    prepared    = ch_main
    stats       = ch_stats

}

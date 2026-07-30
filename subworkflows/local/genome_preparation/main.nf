nextflow.enable.types = true

include { CHECK_GENOME              } from '../../../modules/local/check/genome'
include { SEQKIT_STATS              } from '../../../modules/local/seqkit/stats'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow GENOME_PREPARATION {

    take:
    ch_input

    main:

    // ----------------------------------------------------------
    // CHECK HEADERS OR WHOLE PROTEIN DB AND RAISE ERROR IF UNHANDLED CHARACTERS
    // ----------------------------------------------------------

    CHECK_GENOME( ch_input )
    ch_input = ch_input.join(CHECK_GENOME.out.checked, by: 'id')
    
    // ----------------------------------------------------------
    // COMPUTE STATS ABOUT GENOME
    // ----------------------------------------------------------

    SEQKIT_STATS ( ch_input )
    ch_stats = SEQKIT_STATS.out.stats

    ch_prepared_genome = ch_input
                            .join(ch_stats, by: 'id')
                            .map {
                                rec ->
                                    def csv = rec.genome_stats.splitCsv( header: true, sep: '\t', limit: 1 ).collect()
                                    rec + record(genome_size: csv[0].sum_len.toInteger())
                            }

    emit:
    prepared_genome     = ch_prepared_genome
    stats               = ch_stats

}

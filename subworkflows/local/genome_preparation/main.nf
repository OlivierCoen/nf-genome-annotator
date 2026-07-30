nextflow.enable.types = true

include { CHECK_GENOME              } from '../../../modules/local/check/genome'
include { SEQKIT_STATS              } from '../../../modules/nf-core/seqkit/stats'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow GENOME_PREPARATION {

    take:
    ch_input

    main:
ch_input.view{ v -> "1 $v"}
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
                            .join( ch_stats )
                            .map {
                                meta, genome, stats ->
                                    def csv = stats.splitCsv( header: true, sep: '\t', limit: 1 ).collect()
                                    def stat_row = csv[0]
                                    [ meta + [ genome_size: stat_row.sum_len ], genome]
                            }

    emit:
    prepared_genome     = channel.empty()
    stats = channel.empty()
    //stats               = ch_stats

}

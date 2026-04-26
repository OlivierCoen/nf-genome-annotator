include { CHECK_FASTA              } from '../../../modules/local/check_fasta'
include { SEQKIT_STATS              } from '../../../modules/nf-core/seqkit/stats'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow GENOME_PREPARATION {

    take:
    ch_genome

    main:

    // ----------------------------------------------------------
    // CHECK HEADERS OR WHOLE PROTEIN DB AND RAISE ERROR IF UNHANDLED CHARACTERS
    // ----------------------------------------------------------

    CHECK_FASTA(
        ch_genome,
        type="dna",
        fix_headers=false,
        fix_sequences=false
    )


    // ----------------------------------------------------------
    // COMPUTE STATS ABOUT GENOME
    // ----------------------------------------------------------

    SEQKIT_STATS ( ch_genome )

    ch_prepared_genome = ch_genome
                            .join( SEQKIT_STATS.out.stats )
                            .map {
                                meta, genome, stats ->
                                    def csv = stats.splitCsv( header: true, sep: '\t', limit: 1 ).collect()
                                    def stat_row = csv[0]
                                    [ meta + [ genome_size: stat_row.sum_len ], genome]
                            }

    emit:
    prepared_genome     = ch_prepared_genome

}

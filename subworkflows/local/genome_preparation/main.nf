nextflow.enable.types = true

include { CHECK_GENOME              } from '../../../modules/local/check/genome'
include { SEQKIT_STATS              } from '../../../modules/local/seqkit/stats'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Genome {
    id: String
    fasta: Path
}


workflow GENOME_PREPARATION {

    take:
    ch_input: Channel<Genome>

    main:

    // ----------------------------------------------------------
    // CHECK HEADERS OR WHOLE PROTEIN DB AND RAISE ERROR IF UNHANDLED CHARACTERS
    // ----------------------------------------------------------

    CHECK_GENOME( ch_input )
    ch_prepared = ch_input.join(CHECK_GENOME.out, by: 'id')
    
    // ----------------------------------------------------------
    // COMPUTE STATS ABOUT GENOME
    // ----------------------------------------------------------

    SEQKIT_STATS( ch_prepared )

    ch_prepared = ch_prepared
                .join(SEQKIT_STATS.out, by: 'id')
                .map {
                    rec ->
                        def csv = rec.stats.splitCsv( header: true, sep: '\t', limit: 1 ).collect()
                        rec + record(genome_size: csv[0].sum_len.toInteger())
                }

    emit:
    prepared = ch_prepared

}

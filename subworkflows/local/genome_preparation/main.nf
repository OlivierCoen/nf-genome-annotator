nextflow.enable.types = true

include { CHECK_GENOME              } from '../../../modules/local/check/genome'
include { SEQKIT_STATS              } from '../../../modules/local/seqkit/stats'
include { GET_TAXONOMY              } from '../../../modules/local/get_taxonomy'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Genome {
    id: String
    fasta: Path
    species: String
}


workflow GENOME_PREPARATION {

    take:
    ch_input: Channel<Genome>

    main:

    // ----------------------------------------------------------
    // CHECK HEADERS OR WHOLE PROTEIN DB AND RAISE ERROR IF UNHANDLED CHARACTERS
    // ----------------------------------------------------------

    ch_checked_genome = CHECK_GENOME( ch_input )
    ch_prepared = ch_input.join(ch_checked_genome, by: 'id')
    
    // ----------------------------------------------------------
    // COMPUTE STATS ABOUT GENOME
    // ----------------------------------------------------------

    ch_stats = SEQKIT_STATS( ch_prepared )

    ch_prepared = ch_prepared
                    .join(ch_stats, by: 'id')
                    .map {
                        rec ->
                            def csv = rec.stats.splitCsv( header: true, sep: '\t', limit: 1 ).collect()
                            rec + record(genome_size: csv[0].sum_len.toInteger())
                    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FETCH NCBI TAXON ID AND LINEAGE
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_taxonomy = GET_TAXONOMY(
        ch_prepared.map{ rec -> rec.species }.unique()
    )

    emit:
    ch_prepared.join( ch_taxonomy, by: 'species' )

}

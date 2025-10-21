include { SEQKIT_STATS                                              } from '../../../modules/nf-core/seqkit/stats'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow GENOME_PREPARATION {

    take:
    ch_genome

    main:

    SEQKIT_STATS ( ch_genome )

    ch_genome
         .join( SEQKIT_STATS.out.stats )
         .map {
            meta, genome, stats ->
                def csv = stats.splitCsv( header: true, sep: '\t', limit: 1 ).collect()
                def stat_row = csv[0]
                [ meta + [ genome_size: stat_row.sum_len ], genome]
        }
        .set { ch_prepared_genome }

    emit:
    prepared_genome     = ch_prepared_genome

}

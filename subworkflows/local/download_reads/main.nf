nextflow.enable.types = true

include { DOWNLOAD_SRA           } from '../download_sra'
include { DOWNLOAD_ENA           } from '../download_ena'

// ----------------------------------------------------------------------------
// DOWNLOAD READS FROM PUBLIC DATABASES
// ----------------------------------------------------------------------------

record ExperimentIDs {
    id: String
    rnaseq_experiment_ids: Iterable<String>
}

workflow DOWNLOAD_READS {

    take:
    ch_input: Channel<ExperimentIDs>

    main:

    // creating a channel containing unique public ids (SRA / ENA)
    // for both short reads and long reads
    ch_all_sra_ids = ch_input
                        .map { rec -> rec.short_read_sra_ids + rec.long_read_sra_ids }
                        .flatMap{ sra_ids -> sra_ids.collect() }
                        .filter{ sra_id -> sra_id != null }
                        .unique()

    // ------------------------------------------------------------------------------------
    // DOWNLOAD SRA DATA
    // ------------------------------------------------------------------------------------

    ch_downloaded_sra = DOWNLOAD_SRA(
        ch_all_sra_ids.filter{ id -> id.startsWith('SR') || id.startsWith('DR') }
    )

    // ------------------------------------------------------------------------------------
    // DOWNLOAD ENA DATA
    // ------------------------------------------------------------------------------------

    ch_downloaded_ena = DOWNLOAD_ENA(
        ch_all_sra_ids.filter{ id -> id.startsWith('ER') }
    )

    // ------------------------------------------------------------------------------------
    // ASSOCIATE DOWNLOADED READS BACK TO SAMPLE IDS
    // ------------------------------------------------------------------------------------

    ch_downloaded_reads = ch_downloaded_sra
                            .mix( ch_downloaded_ena )
                            .map{ rec -> record(sra_id: rec.id, reads: rec.reads) }
                        
    // associating back to the corresponding sample IDs
    // TODO: simplify when groupBy can handle records

    // SHORT READS
    ch_downloaded_short_reads = ch_input
                                    .flatMap{
                                        rec -> rec.short_read_sra_ids.collect{ value -> record(id: rec.id, sra_id: value) }
                                    }
                                    .join(ch_downloaded_reads, by: 'sra_id')
                                    .map{ rec -> tuple(rec.id, [rec.sra_id, rec.reads]) }
                                    .groupTuple()
                                    .map{ id, tuples ->
                                        record(
                                            id: id,
                                            downloaded_short_reads: tuples.collect{ tup ->
                                                record(id: tup[0], reads: tup[1].flatten())
                                            }
                                        )
                                    }

    // LONG READS
    ch_downloaded_long_reads = ch_input
                                    .flatMap{
                                        rec -> rec.long_read_sra_ids.collect{ value -> record(id: rec.id, sra_id: value) }
                                    }
                                    .join(ch_downloaded_reads, by: 'sra_id')
                                    .map{ rec -> tuple(rec.id, [rec.sra_id, rec.reads]) }
                                    .groupTuple()
                                    .map{ id, tuples ->
                                        record(
                                            id: id,
                                            downloaded_long_reads: tuples.collect{ tup ->
                                                record(id: tup[0], reads: tup[1].flatten())
                                            }
                                        )
                                    }

    ch_input = ch_input
                .join(ch_downloaded_short_reads, by: 'id', remainder: true)
                .join(ch_downloaded_long_reads, by: 'id', remainder: true)

    emit:
    ch_input

}

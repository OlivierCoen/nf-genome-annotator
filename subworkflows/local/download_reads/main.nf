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
    ch_ids: Channel<ExperimentIDs>

    main:

    // creating a channel containing unique public ids (SRA / ENA)
    ch_experiment_ids = ch_ids.flatMap{ rec -> rec.rnaseq_experiment_ids.collect() }.unique()

    // ------------------------------------------------------------------------------------
    // DOWNLOAD SRA DATA
    // ------------------------------------------------------------------------------------

    DOWNLOAD_SRA(
        ch_experiment_ids.filter{ id -> id.startsWith('SR') || id.startsWith('DR') }
    )

    // ------------------------------------------------------------------------------------
    // DOWNLOAD ENA DATA
    // ------------------------------------------------------------------------------------

    DOWNLOAD_ENA(
        ch_experiment_ids.filter{ id -> id.startsWith('ER') }
    )

    // ------------------------------------------------------------------------------------
    // ASSOCIATE DOWNLOADED READS BACK TO SAMPLE IDS
    // ------------------------------------------------------------------------------------

    ch_downloaded_reads = DOWNLOAD_SRA.out.reads
                            .mix( DOWNLOAD_ENA.out.reads )
                            .map{ rec -> record(experiment_id: rec.id, reads: rec.reads) }

    // associating back to the corresponding sample IDs
    // TODO: simplify when groupBy can handle records
    ch_reads = ch_ids
                .flatMap{
                    rec -> rec.rnaseq_experiment_ids.collect{ value -> record(id: rec.id, experiment_id: value) }
                }
                .join(ch_downloaded_reads, by: 'experiment_id')
                .map{ rec -> tuple(rec.id, [rec.experiment_id, rec.reads]) }
                .groupTuple()
                .map{ id, tuples ->
                    record(
                        id: id,
                        downloaded_rnaseq_fastqs: tuples.collect{ tup ->
                            record(id: tup[0], reads: tup[1].flatten())
                        }
                    )
                }

    emit:
    reads = ch_reads

}

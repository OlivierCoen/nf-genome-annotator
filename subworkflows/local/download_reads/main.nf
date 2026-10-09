nextflow.enable.types = true

include { DOWNLOAD_SRA                              } from '../download_sra'
include { DOWNLOAD_ENA                              } from '../download_ena'

include { SEQTK_SAMPLE as SAMPLE_READS              } from '../../../modules/local/seqtk/sample'

// ----------------------------------------------------------------------------
// DOWNLOAD READS FROM PUBLIC DATABASES
// ----------------------------------------------------------------------------

record ExperimentIDs {
    id: String
    rnaseq_experiment_ids: List<String>
}

workflow DOWNLOAD_READS {

    take:
    ch_input: Channel<ExperimentIDs>
    read_sampling_size: Float

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
    // SAMPLE READS
    // ------------------------------------------------------------------------------------

    ch_downloaded_reads = ch_downloaded_sra.mix( ch_downloaded_ena )

    if ( read_sampling_size ) {

        ch_sampled_reads = SAMPLE_READS(
            ch_downloaded_reads,
            read_sampling_size
        )
        ch_downloaded_reads = ch_downloaded_reads.join( ch_sampled_reads, by: 'id' )
    
    }
    
    // ------------------------------------------------------------------------------------
    // ASSOCIATE DOWNLOADED READS BACK TO SAMPLE IDS
    // ------------------------------------------------------------------------------------
     
    // associating back to the corresponding sample IDs
    // TODO: simplify when groupBy can handle records

    ch_downloaded_reads = ch_downloaded_reads.map{ rec -> record(sra_id: rec.id, reads: rec.reads) }

    // SHORT READS
    ch_downloaded_short_reads = ch_input
                                    .flatMap{ // getting all SRA Ids to prepare the join
                                        rec -> rec.short_read_sra_ids.collect{ value -> record(id: rec.id, sra_id: value) }
                                    }
                                    .join(ch_downloaded_reads, by: 'sra_id')
                                    .map{ rec -> tuple(rec.id, record(id: rec.sra_id, reads: rec.reads)) }
                                    .groupBy()
                                    .map{ id, records ->  record(id: id, downloaded_short_reads: records) }

    // LONG READS
    ch_downloaded_long_reads = ch_input
                                    .flatMap{ // getting all SRA Ids to prepare the join
                                        rec -> rec.long_read_sra_ids.collect{ value -> record(id: rec.id, sra_id: value) }
                                    }
                                    .join(ch_downloaded_reads, by: 'sra_id')
                                    .map{ rec -> tuple(rec.id, record(id: rec.sra_id, fastq: rec.reads[0])) } // changing key: 'reads' -> 'fastq'
                                    .groupBy()
                                    .map{ id, records ->  record(id: id, downloaded_long_reads: records) }
                                    

    ch_input = ch_input
                .join(ch_downloaded_short_reads, by: 'id', remainder: true)
                .join(ch_downloaded_long_reads, by: 'id', remainder: true)

    emit:
    ch_input

}

nextflow.enable.types = true

include { FASTPLONG                             } from '../../../modules/local/fastplong'
include { FASTQC as FASTQC_RAW                  } from '../../../modules/local/fastqc'
include { FASTQC as FASTQC_CLEANED              } from '../../../modules/local/fastqc'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Read {
    id: String
    reads: List<Path>
}

record Input {
    id: String
    long_reads: List<Read>
}


workflow LONG_READ_PREPARATION {

    take:
    ch_input: Channel
    skip_fastqc: Boolean
    skip_fastqc_raw: Boolean
    skip_fastqc_cleaned: Boolean
    skip_long_read_cleaning: Boolean

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // EXTRACT READS 
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // make channel containing one pair of fastq files (or one single fastq file in case of single end)
    // per element
    ch_reads = ch_input
                .flatMap{ rec -> rec.long_reads.collect{ subrec ->
                    record(
                        id: rec.id,
                        read_id: subrec.id,
                        reads: subrec.reads
                    ) }
                }

    // ---------------------------------------------------------------------
    // QUALITY CONTROL ON RAW READS
    // ---------------------------------------------------------------------

    if ( !skip_fastqc && !skip_fastqc_raw ) {
        FASTQC_RAW ( ch_reads )
    }

    // ---------------------------------------------------------------------
    // CLEANING (TRIMMING + ADAPTER REMOVAL) + QC ON CLEANED READS
    // ---------------------------------------------------------------------

    if ( !skip_long_read_cleaning ) {

        ch_cleaned_reads = FASTPLONG(
            ch_reads.map { rec -> record(id: rec.id, fastq: rec.reads[0]) }
        )
        ch_reads = ch_reads.join( ch_cleaned_reads, by: 'id' )
    
            if ( !skip_fastqc && !skip_fastqc_cleaned ) {
                FASTQC_CLEANED ( ch_reads )
            }

    }

    // NOTE: for now, no structural annotator integrated in the pipeline accepts Isoseq data in BAM format
    // if Braker4 gets integrated in the pipeline in the future, we'll need to add mapping steps here

    // ---------------------------------------------------------------------
    // FORMAT OUTPUT CHANNEL
    // ---------------------------------------------------------------------

    ch_reads = ch_reads
                .map { rec -> tuple(rec.sample_id, rec.fastq) }
                .groupTuple()
                .map { id, fastqs -> record(id: id, new_long_reads: fastqs) }

    emit:
    ch_reads
}

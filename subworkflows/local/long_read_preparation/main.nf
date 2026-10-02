nextflow.enable.types = true

include { FASTPLONG                             } from '../../../modules/local/fastplong'
include { FASTQC as FASTQC_RAW                  } from '../../../modules/local/fastqc'
include { FASTQC as FASTQC_CLEANED              } from '../../../modules/local/fastqc'

include { FASTQ_ALIGN_MINIMAP2                  } from '../fastq_align_minimap2'


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
            ch_reads.map { rec -> record(id: rec.id, read_id: rec.read_id, fastq: rec.reads[0]) }
        )
        ch_reads = ch_reads.join( ch_cleaned_reads, by: 'read_id' )
    
            if ( !skip_fastqc && !skip_fastqc_cleaned ) {
                FASTQC_CLEANED ( ch_reads )
            }

    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP READS TO GENOME
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_mapped = FASTQ_ALIGN_MINIMAP2( ch_reads )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GROUP READS & BAM FILES BY SAMPLE ID
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_reads = ch_reads
                .map { rec -> tuple(rec.id, rec.read_id, rec.fastq) }
                .groupTuple()
                .map { id, read_ids, fastqs -> 
                    def new_long_reads = [read_ids, fastqs].collect { read_id, fastq -> 
                        record(read_id: read_id, fastq: fastqs) 
                    }
                    record(id: id, new_long_reads: new_long_reads)
                }

    ch_mapped = ch_mapped
                .map { rec -> tuple(rec.id, rec.read_id, rec.bam) }
                .groupTuple()
                .map { id, read_ids, bams -> 
                    def new_long_read_bams = [read_ids, bams].collect { read_id, bam -> 
                        record(id: read_id, bam: bam) 
                    }
                    record(id: id, new_long_read_bams: new_long_read_bams) 
                }

    emit:
    ch_reads.join( ch_mapped, by: 'id' )

}

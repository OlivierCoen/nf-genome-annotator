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
    read_id: String
    fastq: List<Path>
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
                        fasta: rec.fasta,
                        read_id: subrec.id,
                        fastq: subrec.fastq
                    ) }
                }

    // ---------------------------------------------------------------------
    // QUALITY CONTROL ON RAW READS
    // ---------------------------------------------------------------------

    if ( !skip_fastqc && !skip_fastqc_raw ) {
        FASTQC_RAW ( 
            ch_reads.map { rec -> record(id: rec.id, read_id: rec.read_id, reads: [rec.fastq]) }
        )
    }

    // ---------------------------------------------------------------------
    // CLEANING (TRIMMING + ADAPTER REMOVAL) + QC ON CLEANED READS
    // ---------------------------------------------------------------------

    if ( !skip_long_read_cleaning ) {

        ch_cleaned_reads = FASTPLONG( ch_reads )
        ch_reads = ch_reads.join( ch_cleaned_reads, by: 'read_id' )
    
            if ( !skip_fastqc && !skip_fastqc_cleaned ) {
                FASTQC_CLEANED ( 
                    ch_reads.map { rec -> record(id: rec.id, read_id: rec.read_id, reads: [rec.fastq]) }
                )
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
                .map { rec -> tuple(rec.id, rec.subMap(['read_id', 'fastq'])) }
                .groupBy()
                .map { id, records -> record(id: id, new_long_reads: records) }

    ch_mapped = ch_mapped
                .map { rec -> tuple(rec.id, rec.subMap(['read_id', 'bam'])) }
                .groupBy()
                .map { id, records -> record(id: id, new_long_read_bams: records) }

    emit:
    ch_reads.join( ch_mapped, by: 'id' )

}

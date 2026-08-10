nextflow.enable.types = true

//
// Read QC, UMI extraction and trimming
//
include { FASTQC as FASTQC_RAW  } from '../../../modules/local/fastqc'
include { FASTQC as FASTQC_TRIM } from '../../../modules/local/fastqc'
include { UMITOOLS_EXTRACT      } from '../../../modules/local/umitools/extract'
include { FASTP                 } from '../../../modules/local/fastp'

record Reads {
    id: String
    reads:  List<Path>
}


workflow FASTQ_FASTQC_UMITOOLS_FASTP {

    take:
    ch_reads: Channel<Reads>
    skip_fastqc: Boolean
    skip_umi_extract: Boolean
    skip_trimming: Boolean

    main:

    if (!skip_fastqc) {
        FASTQC_RAW( ch_reads )
    }

    if ( !skip_umi_extract ) {
        UMITOOLS_EXTRACT( ch_reads )
        ch_reads = UMITOOLS_EXTRACT.out
    }

    if ( !skip_trimming ) {

        FASTP( ch_reads )
        ch_reads = ch_reads.join( FASTP.out, by: 'id' )

        if ( !skip_fastqc ) {
            FASTQC_TRIM( ch_reads )
        }
    }

    emit:
    reads = ch_reads
}

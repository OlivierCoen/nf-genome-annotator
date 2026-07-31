nextflow.enable.types = true

include { FETCH_ENA_FASTQ_URLS        } from '../../../modules/local/ena/fetch_ena_fastq_urls'
include { DOWNLOAD_ENA_FASTQ          } from '../../../modules/local/ena/download_ena_fastq'

// ----------------------------------------------------------------------------
// DOWNLOAD FASTQ SEQUENCING READS FROM THE ENA ARCHIVE.
// ----------------------------------------------------------------------------

workflow DOWNLOAD_ENA {

    take:
    ch_ids: Channel<String>

    main:

    // ----------------------------------------
    // FETCH ENA FASTQ URLS
    // ----------------------------------------

    FETCH_ENA_FASTQ_URLS ( ch_ids )

    // ---------------------------------------------------------------
    // DOWNLOAD ENA FASTQ FILES
    // ---------------------------------------------------------------

    DOWNLOAD_ENA_FASTQ ( FETCH_ENA_FASTQ_URLS.out )


    emit:
    reads = DOWNLOAD_ENA_FASTQ.out

}

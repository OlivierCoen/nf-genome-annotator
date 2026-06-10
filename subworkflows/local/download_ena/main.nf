include { FETCH_ENA_FASTQ_URLS        } from '../../../modules/local/ena/fetch_ena_fastq_urls'
include { DOWNLOAD_ENA_FASTQ          } from '../../../modules/local/ena/download_ena_fastq'

// ----------------------------------------------------------------------------
// DOWNLOAD FASTQ SEQUENCING READS FROM THE ENA ARCHIVE.
// ----------------------------------------------------------------------------

workflow DOWNLOAD_ENA {
    take:
    ch_ena_ids   // channel: [ val(meta), val(ena_id) ]

    main:

    // ----------------------------------------
    // FETCH ENA FASTQ URLS
    // ----------------------------------------

    FETCH_ENA_FASTQ_URLS ( ch_ena_ids )

    // ---------------------------------------------------------------
    // DOWNLOAD ENA FASTQ FILES
    // ---------------------------------------------------------------

    DOWNLOAD_ENA_FASTQ ( FETCH_ENA_FASTQ_URLS.out.ftp_urls)


    emit:
    reads               = DOWNLOAD_ENA_FASTQ.out.fastq

}

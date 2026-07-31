nextflow.enable.types = true

include { CUSTOM_SRATOOLSNCBISETTINGS } from '../../../modules/nf-core/custom/sratoolsncbisettings'
include { SRATOOLS_PREFETCH           } from '../../../modules/local/sratools/prefetch'
include { SRATOOLS_FASTERQDUMP        } from '../../../modules/local/sratools/fasterqdump'

// ----------------------------------------------------------------------------
// DOWNLOAD FASTQ SEQUENCING READS FROM THE NCBI'S SEQUENCE READ ARCHIVE (SRA).
// ----------------------------------------------------------------------------

workflow DOWNLOAD_SRA {

    take:
    ch_ids: Channel<String>

    main:

    // --------------------------------------------------------
    // DETECT EXISTING NCBI USER SETTINGS OR CREATE NEW ONES.
    // --------------------------------------------------------

    CUSTOM_SRATOOLSNCBISETTINGS( [] )
    ch_ncbi_settings = CUSTOM_SRATOOLSNCBISETTINGS.out.ncbi_settings

    // ----------------------------------------
    // PREFETCH SEQUENCING READS IN SRA FORMAT.
    // ----------------------------------------

    SRATOOLS_PREFETCH (
        ch_ids,
        ch_ncbi_settings
    )

    ch_sra = SRATOOLS_PREFETCH.out
                .flatMap { rec -> // transpose (typed version): when multiple SRRs are downloaded for a specific SRA ID, we split them
                    rec.sra.collect { value -> record(id: rec.id, sra: value) }
                }

    // ---------------------------------------------------------------
    // CONVERT THE SRA FORMAT INTO ONE OR MORE COMPRESSED FASTQ FILES.
    // ---------------------------------------------------------------

    SRATOOLS_FASTERQDUMP (
        ch_sra,
        ch_ncbi_settings
    )

    emit:
    reads = SRATOOLS_FASTERQDUMP.out
}

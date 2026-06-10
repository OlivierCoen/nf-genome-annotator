include { CUSTOM_SRATOOLSNCBISETTINGS } from '../../../modules/nf-core/custom/sratoolsncbisettings'
include { SRATOOLS_PREFETCH           } from '../../../modules/local/sratools/prefetch'
include { SRATOOLS_FASTERQDUMP        } from '../../../modules/local/sratools/fasterqdump'

// ----------------------------------------------------------------------------
// DOWNLOAD FASTQ SEQUENCING READS FROM THE NCBI'S SEQUENCE READ ARCHIVE (SRA).
// ----------------------------------------------------------------------------

workflow DOWNLOAD_SRA {
    take:
    ch_sra_ids   // channel: [ val(meta), val(sra_id) ]

    main:

    // --------------------------------------------------------
    // DETECT EXISTING NCBI USER SETTINGS OR CREATE NEW ONES.
    // --------------------------------------------------------

    CUSTOM_SRATOOLSNCBISETTINGS ( ch_sra_ids.collect() )
    ch_ncbi_settings = CUSTOM_SRATOOLSNCBISETTINGS.out.ncbi_settings

    // ----------------------------------------
    // PREFETCH SEQUENCING READS IN SRA FORMAT.
    // ----------------------------------------

    SRATOOLS_PREFETCH (
        ch_sra_ids,
        ch_ncbi_settings
    )

    ch_sra = SRATOOLS_PREFETCH.out.sra
                .transpose() // when multiple SRRs are downloaded for a specific SRA ID, we split them
                .map {
                    meta, sra ->
                        def new_meta = [ id: sra.name ] + meta
                        [ new_meta, sra ]
                }

    // ---------------------------------------------------------------
    // CONVERT THE SRA FORMAT INTO ONE OR MORE COMPRESSED FASTQ FILES.
    // ---------------------------------------------------------------

    SRATOOLS_FASTERQDUMP (
        ch_sra,
        ch_ncbi_settings
    )
    ch_sra_reads = SRATOOLS_FASTERQDUMP.out.reads


    emit:
    reads               = ch_sra_reads
}

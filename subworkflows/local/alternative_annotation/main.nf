include { AGAT_SPKEEPLONGESTISOFORM as AGAT_KEEP_LONGEST_ISOFORM } from '../../../modules/local/agat/spkeeplongestisoform'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow ALTERNATIVE_ANNOTATIONS {

    take:
    ch_gff

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // KEEPING ONLY LONGEST ISOFORMS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    AGAT_KEEP_LONGEST_ISOFORM ( ch_gff )

    emit:
    longest_isoforms_gff = AGAT_KEEP_LONGEST_ISOFORM.out.gff

}

include { REPEATMODELER_BUILDDATABASE as BUILDDATABASE              } from '../../../modules/nf-core/repeatmodeler/builddatabase'
include { REPEATMODELER_REPEATMODELER as REPEATMODELER              } from '../../../modules/nf-core/repeatmodeler/repeatmodeler'
include { REPEATMASKER_REPEATMASKER as REPEATMASKER                 } from '../../../modules/local/repeatmasker/repeatmasker'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow SOFTMASKING {

    take:
    ch_genome

    main:

    ch_versions = Channel.empty()

    BUILDDATABASE( ch_genome )

    REPEATMODELER( BUILDDATABASE.out.db )

    ch_genome
        .join ( REPEATMODELER.out.fasta )
        .set { ch_repeat_masker_input }

    REPEATMASKER( ch_repeat_masker_input )

    ch_versions
        .mix ( BUILDDATABASE.out.versions )
        .mix ( REPEATMODELER.out.versions )
        .set { ch_versions }


    emit:
    softmasked_genome       = REPEATMASKER.out.masked
    versions                = ch_versions

}


include { AGAT_SQREMOVEREDUNDANTENTRIES                } from '../../../modules/local/agat/sq_remove_redundant_entries'
include { AGAT_SPFIXOVERLAPINGGENES                    } from '../../../modules/local/agat/spfixoverlapinggenes'
include { AGAT_SPKEEPLONGESTISOFORM                    } from '../../../modules/nf-core/agat/spkeeplongestisoform'
include { AGAT_SPFILTERINCOMPLETEGENECODINGMODELS      } from '../../../modules/local/agat/spfilterincompletegenecodingmodels'
// include { AGAT_SPFIXLONGESTORF                         } from '../../../modules/local/agat/spfixlongestorf'
include { AGAT_SPFIXCDSPHASES                          } from '../../../modules/local/agat/spfixcdsphases'
include { AGAT_SPMANAGEIDS                             } from '../../../modules/local/agat/spmanageids'
include { AGAT_CONVERTSPGFF2GTF                        } from '../../../modules/nf-core/agat/convertspgff2gtf'
include { AGAT_SPSTATISTICS                            } from '../../../modules/nf-core/agat/spstatistics'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow CLEAN_ANNOTATIONS {

    take:
    ch_gtf
    ch_genome

    main:

    ch_versions = Channel.empty()

    AGAT_SQREMOVEREDUNDANTENTRIES ( ch_gtf, [] )

    AGAT_SPKEEPLONGESTISOFORM ( AGAT_SQREMOVEREDUNDANTENTRIES.out.gff, [] )

    AGAT_SPFIXOVERLAPINGGENES ( AGAT_SPKEEPLONGESTISOFORM.out.gff, [] )

    AGAT_SPFILTERINCOMPLETEGENECODINGMODELS ( AGAT_SPFIXOVERLAPINGGENES.out.gff.join( ch_genome ), [] )

    AGAT_SPFIXCDSPHASES ( AGAT_SPFILTERINCOMPLETEGENECODINGMODELS.out.gff.join( ch_genome ), [] )

    AGAT_SPMANAGEIDS ( AGAT_SPFIXCDSPHASES.out.gff, [] )

    AGAT_CONVERTSPGFF2GTF ( AGAT_SPMANAGEIDS.out.gff )
    AGAT_CONVERTSPGFF2GTF.out.output_gtf.set { ch_cleaned_gtf }

    AGAT_SPSTATISTICS ( ch_cleaned_gtf )


    ch_versions
        .mix ( AGAT_CONVERTSPGFF2GTF.out.versions )
        .mix ( AGAT_SPKEEPLONGESTISOFORM.out.versions )
        .mix ( AGAT_SPSTATISTICS.out.versions )
        .set { ch_versions }


    emit:
    gtf                     = ch_cleaned_gtf
    stats                   = AGAT_SPSTATISTICS.out.stats_txt
    versions                = ch_versions

}


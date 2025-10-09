include { AGAT_SQREMOVEREDUNDANTENTRIES                } from '../../../modules/local/agat/sq_remove_redundant_entries'
include { AGAT_SPFIXOVERLAPINGGENES                    } from '../../../modules/local/agat/spfixoverlapinggenes'
include { AGAT_SPKEEPLONGESTISOFORM                    } from '../../../modules/local/agat/spkeeplongestisoform'
include { AGAT_SPFILTERINCOMPLETEGENECODINGMODELS      } from '../../../modules/local/agat/spfilterincompletegenecodingmodels'
// include { AGAT_SPFIXLONGESTORF                         } from '../../../modules/local/agat/spfixlongestorf'
include { AGAT_SPFIXCDSPHASES                          } from '../../../modules/local/agat/spfixcdsphases'
include { AGAT_SPMANAGEIDS                             } from '../../../modules/local/agat/spmanageids'
include { AGAT_CONVERTSPGFF2GTF                        } from '../../../modules/local/agat/convertspgff2gtf'



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

    AGAT_SQREMOVEREDUNDANTENTRIES ( ch_gtf, [] )

    AGAT_SPKEEPLONGESTISOFORM ( AGAT_SQREMOVEREDUNDANTENTRIES.out.gff, [] )

    AGAT_SPFIXOVERLAPINGGENES ( AGAT_SPKEEPLONGESTISOFORM.out.gff, [] )

    AGAT_SPFILTERINCOMPLETEGENECODINGMODELS ( AGAT_SPFIXOVERLAPINGGENES.out.gff.join( ch_genome ), [] )

    AGAT_SPFIXCDSPHASES ( AGAT_SPFILTERINCOMPLETEGENECODINGMODELS.out.gff.join( ch_genome ), [] )

    AGAT_SPMANAGEIDS ( AGAT_SPFIXCDSPHASES.out.gff, [] )

    AGAT_CONVERTSPGFF2GTF ( AGAT_SPMANAGEIDS.out.gff )


    emit:
    gtf                     = AGAT_CONVERTSPGFF2GTF.out.output_gtf

}


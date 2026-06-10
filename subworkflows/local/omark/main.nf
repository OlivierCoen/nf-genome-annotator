include { OMARK_DOWNLOADDB                       } from '../../../modules/local/omark/download_db'
include { OMARK_OMAMERSEARCH                     } from '../../../modules/local/omark/omamer_search'
include { OMARK_EXTRACT_TRANSCRIPT_ISOFORMS      } from '../../../modules/local/omark/extract_transcript_isoforms'
include { OMARK_OMARK                            } from '../../../modules/local/omark/omark'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow OMARK {

    take:
    ch_proteome
    ch_gff
    omamer_db_url

    main:

    OMARK_DOWNLOADDB( omamer_db_url )
    ch_omark_db = OMARK_DOWNLOADDB.out.db.collect()
    
    OMARK_OMAMERSEARCH(
        ch_proteome,
        ch_omark_db
    )

    OMARK_EXTRACT_TRANSCRIPT_ISOFORMS( ch_gff )

    ch_omark_input = ch_proteome
                       .join( OMARK_OMAMERSEARCH.out.omamer ) 
                       .join( OMARK_EXTRACT_TRANSCRIPT_ISOFORMS.out.isoforms )

    OMARK_OMARK( 
        ch_omark_input,
        ch_omark_db
    )
}

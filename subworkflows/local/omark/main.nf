nextflow.enable.types = true

include { OMARK_DOWNLOADDB                       } from '../../../modules/local/omark/download_db'
include { OMARK_OMAMERSEARCH                     } from '../../../modules/local/omark/omamer_search'
include { OMARK_EXTRACT_TRANSCRIPT_ISOFORMS      } from '../../../modules/local/omark/extract_transcript_isoforms'
include { OMARK_OMARK                            } from '../../../modules/local/omark/omark'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    gff: Path
    busco_lineage: String
    proteome: Path
    other_proteomes: Iterable<Path>
}

workflow OMARK {

    take:
    ch_input: Channel<Input>
    omamer_db_url: String
    omamer_db: String

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DOWNLOAD DB
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( omamer_db ) {
        //log.info "Using the provided OMAMER database: ${omamer_db}"
        ch_omark_db = channel.fromPath( omamer_db, checkExists: true )
    } else {
        ch_omark_db = OMARK_DOWNLOADDB( omamer_db_url )
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PREPARE OMARK RUN
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ch_omamer_search_out = OMARK_OMAMERSEARCH(
        ch_input.map { rec -> record(id: rec.id, fasta: rec.proteome) },
        ch_omark_db
    )

    ch_transcript_isoforms = OMARK_EXTRACT_TRANSCRIPT_ISOFORMS( ch_input )

    ch_input = ch_input
                .join( ch_omamer_search_out, by: 'id' )
                .join( ch_transcript_isoforms, by: 'id' )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // OMARK
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_omark_out = OMARK_OMARK( 
        ch_input.map { rec -> record(id: rec.id, fasta: rec.proteome, omamer: rec.omamer, transcript_isoforms: rec.transcript_isoforms) },
        ch_omark_db
    )

}

include { DOWNLOAD_SRA           } from '../download_sra'
include { DOWNLOAD_ENA           } from '../download_ena'


// ----------------------------------------------------------------------------
// DOWNLOAD READS FROM PUBLIC DATABASES
// ----------------------------------------------------------------------------

workflow DOWNLOAD_READS {

    take:
    ch_ids   // channel: [ val(meta), val(id) ]

    main:

    ch_db_specific_ids = ch_ids
                            .branch { meta, id ->
                                        sra: id.startsWith('SR')
                                        ena: id.startsWith('ER')
                            }

    // ------------------------------------------------------------------------------------
    // DOWNLOAD SRA DATA
    // ------------------------------------------------------------------------------------
    
    DOWNLOAD_SRA( ch_db_specific_ids.sra )
    
    // ------------------------------------------------------------------------------------
    // DOWNLOAD ENA DATA
    // ------------------------------------------------------------------------------------
    
    DOWNLOAD_ENA( ch_db_specific_ids.ena )
    
    ch_downloaded_reads = DOWNLOAD_SRA.out.reads
                            .mix ( DOWNLOAD_ENA.out.reads )
                            .map{ meta, reads ->
                                def single_end_state = reads instanceof Path ? true : false
                                [ meta + [ single_end: single_end_state ], reads ]
                            }
    
    
    emit:
    reads = ch_downloaded_reads

}
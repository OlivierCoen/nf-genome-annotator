nextflow.enable.types = true

include { DOWNLOAD_SRA           } from '../download_sra'
include { DOWNLOAD_ENA           } from '../download_ena'

// ----------------------------------------------------------------------------
// DOWNLOAD READS FROM PUBLIC DATABASES
// ----------------------------------------------------------------------------

record PublicIDs {
    id: String
    public_id: Iterable<String>
}

workflow DOWNLOAD_READS {

    take:
    ch_ids: Channel<PublicIDs>

    main:

    // creating a channel containing unique public ids (SRA / ENA)
    ch_public_ids = ch_ids.flatMap{ rec -> rec.rnaseq_public_ids.collect() }.unique()
   
    // ------------------------------------------------------------------------------------
    // DOWNLOAD SRA DATA
    // ------------------------------------------------------------------------------------
    
    DOWNLOAD_SRA( 
        ch_public_ids.filter{ id -> id.startsWith('SR') || id.startsWith('DR') }
    )
    
    // ------------------------------------------------------------------------------------
    // DOWNLOAD ENA DATA
    // ------------------------------------------------------------------------------------
    
    DOWNLOAD_ENA( 
        ch_public_ids.filter{ id -> id.startsWith('ER') }
    )

    // ------------------------------------------------------------------------------------
    // ASSOCIATE DOWNLOADED READS BACK TO SAMPLE IDS
    // ------------------------------------------------------------------------------------

    ch_downloaded_reads = DOWNLOAD_SRA.out.reads
                            .mix( DOWNLOAD_ENA.out.reads )
                            .map{ rec -> record(public_id: rec.id, reads: rec.reads) }

    // the Nextflow syntax is not appropriate here...
    ch_reads = ch_ids
                .flatMap{ 
                    rec -> rec.rnaseq_public_ids.collect{ value -> record(id: rec.id, public_id: value) } 
                }
                .join(ch_downloaded_reads, by: 'public_id')
                .map{ rec -> tuple(rec.id, [rec.public_id, rec.reads]) } 
                .groupTuple()
                .map{ id, tuples -> 
                    println tuples
                    
                }
                //.view{ v -> "after $v"}
                            
    // single_end: rec.reads instanceof Path ? true : false 
    // associating back to the corresponding sample IDs
    emit:
    reads = ch_ids

}
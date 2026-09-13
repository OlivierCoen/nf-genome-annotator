nextflow.enable.types = true

include { FASTPLONG                             } from '../../../modules/local/fastplong'
include { FASTQC as FASTQC_RAW                  } from '../../../modules/local/fastqc'
include { FASTQC as FASTQC_CLEANED              } from '../../../modules/local/fastqc'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Reads {
    id: String
    reads:  List<Path>
}


workflow LONG_READ_PREPARATION {

    take:
    ch_input: Channel
    skip_fastqc: Boolean
    skip_fastqc_raw: Boolean
    skip_fastqc_cleaned: Boolean
    skip_long_read_cleaning: Boolean

    main:

    // ---------------------------------------------------------------------
    // QUALITY CONTROL ON RAW READS
    // ---------------------------------------------------------------------

    if ( !skip_fastqc && !skip_fastqc_raw ) {
        FASTQC_RAW ( ch_input )
    }

    // ---------------------------------------------------------------------
    // CLEANING (TRIMMING + ADAPTER REMOVAL) + QC ON CLEANED READS
    // ---------------------------------------------------------------------

    if ( !skip_long_read_cleaning ) {

        ch_cleaned_reads = FASTPLONG( ch_input )
        ch_input = ch_input.join( ch_cleaned_reads, by: 'id' )
    
            if ( !skip_fastqc && !skip_fastqc_cleaned ) {
                FASTQC_CLEANED ( ch_input )
            }

    }

    // NOTE: for now, no structural annotator integrated in the pipeline accepts Isoseq data in BAM format
    // if Braker4 gets integrated in the pipeline in the future, we'll need to add mapping steps here

    emit:
    ch_input
}

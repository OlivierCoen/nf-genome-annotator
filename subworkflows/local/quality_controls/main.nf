nextflow.enable.types = true


include { AGAT_SPSTATISTICS as AGAT_GTF_STATISTICS                    } from '../../../modules/local/agat/spstatistics'
include { AGAT_SPFUNCTIONALSTATISTICS as AGAT_FUNCTIONAL_STATISTICS   } from '../../../modules/local/agat/spfunctionalstatistics'

include { BUSCO                                                       } from '../busco'
include { OMARK                                                       } from '../omark'

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
    other_proteomes: Set<Path>
}

workflow QUALITY_CONTROLS {

    take:
    ch_input: Channel<Input>
    skip_busco: Boolean
    skip_omark: Boolean
    omamer_db_url: Boolean
    omamer_db: Boolean

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // RUN BUSCO ON THE GENOME AND THE DIFFERENT VERSIONS OF THE PREDICTED PROTEOME
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !skip_busco ) {

        BUSCO( ch_input )

    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ASSESSMENT OF ANNOTATION QUALITY WITH OMARK
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !skip_omark ) {

        OMARK(
            ch_input,
            omamer_db_url,
            omamer_db
        )

    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // METRICS OF STRUCTURAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    AGAT_GTF_STATISTICS ( ch_input )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // METRICS OF FUNCTIONAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    AGAT_FUNCTIONAL_STATISTICS( ch_input )

    emit:
    ch_input

}

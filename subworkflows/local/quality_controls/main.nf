nextflow.enable.types = true

include { BUSCO_DOWNLOAD                                              } from '../../../modules/local/busco/download'
include { BUSCO_BUSCO as BUSCO_GENOME                                 } from '../../../modules/local/busco/busco'
include { BUSCO_BUSCO as BUSCO_PROTEOME                               } from '../../../modules/local/busco/busco'
include { AGAT_SPSTATISTICS as AGAT_GTF_STATISTICS                    } from '../../../modules/local/agat/spstatistics'
include { AGAT_SPFUNCTIONALSTATISTICS as AGAT_FUNCTIONAL_STATISTICS   } from '../../../modules/local/agat/spfunctionalstatistics'

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
    other_proteomes: Iterable<Path>
}

workflow QUALITY_CONTROLS {

    take:
    ch_input: Channel<Input>
    skip_omark: Boolean
    omamer_db_url: Boolean
    omamer_db: Boolean

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DOWNLOAD NECESSARY BUSCO DATASETS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_busco_input = ch_input.filter { rec -> rec.busco_lineage != null }

    ch_busco_downloads = BUSCO_DOWNLOAD(
        ch_busco_input.map { rec -> rec.busco_lineage }.unique()
    )

    ch_busco_input = ch_busco_input.join( ch_busco_downloads, by: 'busco_lineage' )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // BUSCO ON GENOME
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
    BUSCO_GENOME (
        ch_busco_input.map { rec -> record(
                id: rec.id, 
                fasta: rec.fasta, 
                lineage: rec.busco_lineage, 
                download_path: rec.busco_download_path
            ) 
        },
        'genome'
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // BUSCO ON ALL PROTEOMES
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_busco_proteome_input = ch_busco_input.flatMap { rec -> 
        def all_proteomes = [rec.proteome] + rec.other_proteomes
        all_proteomes.collect { proteome -> record( 
            id: rec.id, 
            fasta: proteome, 
            lineage: rec.busco_lineage, 
            download_path: rec.busco_download_path
        ) }
        
    }

    BUSCO_PROTEOME (
        ch_busco_proteome_input,
        'proteins'
    )

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

    // AGAT_GTF_STATISTICS ( ch_input )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // METRICS OF FUNCTIONAL ANNOTATION
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //AGAT_FUNCTIONAL_STATISTICS( ch_input )

}

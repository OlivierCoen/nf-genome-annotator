nextflow.enable.types = true

include { BUSCO_LISTDATASETS as LIST_DATASETS                         } from '../../../modules/local/busco/list_datasets'
include { BUSCO_GET_LINEAGE as GET_LINEAGE                            } from '../../../modules/local/busco/get_lineage'
include { BUSCO_DOWNLOAD as DOWNLOAD_DATASET                          } from '../../../modules/local/busco/download'
include { BUSCO_BUSCO as BUSCO_GENOME                                 } from '../../../modules/local/busco/busco'
include { BUSCO_BUSCO as BUSCO_PROTEOME                               } from '../../../modules/local/busco/busco'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    species: String
    proteome: Path
    other_proteomes: Set<Path>
}

workflow BUSCO {

    take:
    ch_input: Channel<Input>
  

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FETCHING LIST OF BUSCO DATASETS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_available_busco_datasets = LIST_DATASETS([])

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SEARCH A MATCHING LINEAGE BETWEEN THE LIST OF AVAILABLE LINEAGES AND THE LINEAGE OF PARENT ORGANISMS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ch_busco_lineage = GET_LINEAGE( 
        ch_input.map { rec -> rec.species }.unique(),
        ch_available_busco_datasets
    )

    ch_input = ch_input.join(ch_busco_lineage, by: 'species')

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DOWNLOAD NECESSARY BUSCO DATASETS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_busco_input = ch_input.filter { rec -> rec.busco_lineage != null }

    ch_busco_downloads = DOWNLOAD_DATASET(
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


}

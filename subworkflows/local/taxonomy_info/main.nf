nextflow.enable.types = true


include { GET_TAXONOMY_INFO                                  } from '../../../modules/local/get_taxonomy_info'
include { ORTHODB_GET_AVAILABLE_CLADES                       } from '../../../modules/local/orthodb/get_available_clades'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow TAXONOMY_INFO {

    take:
    ch_species: Channel<String>

    main:

   

    

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FETCHING TAXONOMY INFO AND AGGREGATING WITH PREVIOUS DATA FETCHED
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GET_TAXONOMY_INFO(
        ch_species,
        BUSCO_LISTDATASETS.out.datasets.collect(),
        ORTHODB_GETCLADES.out.clades.collect()
    )

    emit:
    taxonomy = GET_TAXONOMY_INFO.out.taxonomy
    
}

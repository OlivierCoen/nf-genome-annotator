include { BUSCO_LISTDATASETS                                 } from '../../../modules/local/busco/list_datasets'
include { GET_TAXONOMY_INFO                                  } from '../../../modules/local/get_taxonomy_info'
include { ORTHODB_GETCLADES                                  } from '../../../modules/local/orthodb/get_clades'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow TAXONOMY_INFO {

    take:
    ch_main

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FETCHING LIST OF BUSCO DATASETS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    BUSCO_LISTDATASETS()

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FETCHING LIST OF ORTHODB CLADES
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ORTHODB_GETCLADES()

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FETCHING TAXONOMY INFO AND AGGREGATING WITH PREVIOUS DATA FETCHED
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GET_TAXONOMY_INFO(
        ch_main.map{ rec -> rec.species }.unique(),
        BUSCO_LISTDATASETS.out.datasets.collect(),
        ORTHODB_GETCLADES.out.clades.collect()
    )

    ch_main = ch_main.join(GET_TAXONOMY_INFO.out.taxonomy, by: 'species')

    emit:
    enriched_with_taxonomy = ch_main
    
}

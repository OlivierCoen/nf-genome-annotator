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
    ch_species

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
        ch_species,
        BUSCO_LISTDATASETS.out.datasets.collect(),
        ORTHODB_GETCLADES.out.clades.collect(),
    )

    emit:
    taxid = GET_TAXONOMY_INFO.out.taxid
    busco_dataset = GET_TAXONOMY_INFO.out.busco_dataset
    orthodb_clade = GET_TAXONOMY_INFO.out.orthodb_clade
    
}

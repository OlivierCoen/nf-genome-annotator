include { BUSCO_LISTDATASETS                                 } from '../../../modules/local/busco/list_datasets'
include { GET_TAXONOMY_INFO                                  } from '../../../modules/local/get_taxonomy_info'
include { ORTHODB_GETCLADES                                  } from '../../../modules/local/orthodb/get_clades'

def formatSpecies(species) {
    def formated_species = species.replaceAll(/\s+/, '_')
    return formated_species
}

def spreadToRelatedgMetas( ch_genome, ch_other ) {
    return ch_genome
                .combine( ch_other )
                .filter{ meta, file1, species, file2 -> formatSpecies(meta.species) == species }
                .map{ meta, file1, species, file2 -> [ meta, file2 ] }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow TAXONOMY_INFO {

    take:
    ch_genome

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

    ch_species  = ch_genome
                    .map{ meta, file -> formatSpecies(meta.species) }
                    .unique()
    
    GET_TAXONOMY_INFO(
        ch_species,
        BUSCO_LISTDATASETS.out.datasets.collect(),
        ORTHODB_GETCLADES.out.clades.collect()
    )

    emit:
    taxid          = spreadToRelatedgMetas( ch_genome, GET_TAXONOMY_INFO.out.taxid )
    busco_lineages = spreadToRelatedgMetas( ch_genome, GET_TAXONOMY_INFO.out.busco_dataset )
    orthodb_clade  = spreadToRelatedgMetas( ch_genome, GET_TAXONOMY_INFO.out.orthodb_clade )
    
}

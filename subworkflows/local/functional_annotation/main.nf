include { EGGNOGMAPPER_DOWNLOADDB           } from '../../../modules/local/eggnogmapper/download_db'
include { EGGNOGMAPPER_EMAPPER              } from '../../../modules/local/eggnogmapper/emapper'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_proteome

    main:

    EGGNOGMAPPER_DOWNLOADDB ( )

    EGGNOGMAPPER_EMAPPER(
        ch_proteome,
        EGGNOGMAPPER_DOWNLOADDB.out.eggnog_data_dir
    )

    emit:
    annotations      = EGGNOGMAPPER_EMAPPER.out.annotations
    orthologs        = EGGNOGMAPPER_EMAPPER.out.orthologs
    hits             = EGGNOGMAPPER_EMAPPER.out.hits


}

/*
workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_proteome

    main:

    ch_versions = Channel.empty()

    if ( params.interproscan_db != null ) {

        interproscan_db = file( params.interproscan_db, checkExists: true )

    } else {

        // DOWNLOADING
        ch_db_url = Channel.value([
            [ id: params.interproscan_db_url.tokenize("/")[-1] - '.tar.gz'],
            params.interproscan_db_url
        ])
        PREPARE_INTERPROSCAN_DB ( ch_db_url )

        interproscan_db = PREPARE_INTERPROSCAN_DB.out.db

        ch_versions = ch_versions.mix( PREPARE_INTERPROSCAN_DB.out.versions )
    }

    INTERPROSCAN( ch_proteome, interproscan_db )

    emit:
    tsv                     = INTERPROSCAN.out.tsv
    versions                = ch_versions

}
*/


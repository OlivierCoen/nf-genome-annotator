include { EGGNOGMAPPER_DOWNLOADDB                      } from '../../../modules/local/eggnogmapper/download_db'
include { EGGNOGMAPPER_EMAPPER as EMAPPER              } from '../../../modules/local/eggnogmapper/emapper'

include { INTERPROSCAN_DOWNLOADDB                      } from '../../../modules/local/interproscan/download_db'
include { INTERPROSCAN_INTERPROSCAN as INTERPROSCAN    } from '../../../modules/local/interproscan/interproscan'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_proteome
    ch_gff

    main:

    ch_versions = Channel.empty()

    if ( params.functional_annotator == "eggnogmapper" ) {

        EGGNOGMAPPER_DOWNLOADDB ( )
        EMAPPER(
            ch_proteome.join( ch_gff ),
            EGGNOGMAPPER_DOWNLOADDB.out.eggnog_data_dir
        )

        ch_versions = ch_versions.mix( EGGNOGMAPPER_DOWNLOADDB.out.versions )

    } else { // interproscan

        if ( params.interproscan_db != null ) {

            interproscan_db = file( params.interproscan_db, checkExists: true )

        } else {

            // DOWNLOADING
            ch_db_url = Channel.value([
                [ id: params.interproscan_db_url.tokenize("/")[-1] - '.tar.gz'],
                params.interproscan_db_url
            ])
            INTERPROSCAN_DOWNLOADDB ( ch_db_url )
            interproscan_db = INTERPROSCAN_DOWNLOADDB.out.db
        }

        INTERPROSCAN( ch_proteome, interproscan_db )

        ch_versions = ch_versions.mix( INTERPROSCAN_DOWNLOADDB.out.versions )
    }


    emit:
    versions         = ch_versions

}

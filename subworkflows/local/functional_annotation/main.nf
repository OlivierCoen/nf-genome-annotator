include { EGGNOGMAPPER_DOWNLOADDB                      } from '../../../modules/local/eggnogmapper/download_db'
include { EGGNOGMAPPER_EMAPPER                         } from '../../../modules/local/eggnogmapper/emapper'

include { INTERPROSCAN_DOWNLOADDB                      } from '../../../modules/local/interproscan5/download_db'
include { INTERPROSCAN_INTERPROSCAN as INTERPROSCAN    } from '../../../modules/local/interproscan5/interproscan'

include { COMPLEMENT_GFF3_WITH_INTERPROSCAN            } from '../../../modules/local/complement_gff3_with_interproscan'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_proteome
    ch_gff
    functional_annotators
    interproscan_db
    interproscan_db_url

    main:

    ch_versions = Channel.empty()

    ch_decorated_gff = ch_gff

    if ( "eggnogmapper" in functional_annotators ) {

        EGGNOGMAPPER_DOWNLOADDB ( )
        ch_eggnog_db = EGGNOGMAPPER_DOWNLOADDB.out.eggnog_data_dir

        EGGNOGMAPPER_EMAPPER(
            ch_proteome.join( ch_gff ),
            ch_eggnog_db
        )

        ch_decorated_gff = EGGNOGMAPPER_EMAPPER.out.decorated_gff

    }

    if ( "interproscan" in functional_annotators ) {

        if ( interproscan_db != null ) {

            interproscan_db = file( interproscan_db, checkExists: true )

        } else {

            // DOWNLOADING
            ch_db_url = Channel.value([
                [ id: interproscan_db_url.tokenize("/")[-1] - '.tar.gz'],
                interproscan_db_url
            ])
            INTERPROSCAN_DOWNLOADDB ( ch_db_url )
            interproscan_db = INTERPROSCAN_DOWNLOADDB.out.db
        }

        INTERPROSCAN( ch_proteome, interproscan_db )

        ch_versions = ch_versions.mix( INTERPROSCAN_DOWNLOADDB.out.versions )

        COMPLEMENT_GFF3_WITH_INTERPROSCAN(
            ch_decorated_gff.join( INTERPROSCAN.out.gff3 )
        )
        ch_decorated_gff = COMPLEMENT_GFF3_WITH_INTERPROSCAN.out.gff3

    }


    emit:
    gff              = ch_decorated_gff
    versions         = ch_versions

}

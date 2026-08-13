nextflow.enable.types = true

include { EGGNOGMAPPER_DOWNLOADDB                      } from '../../../modules/local/eggnogmapper/download_db'
include { EGGNOGMAPPER_EMAPPER                         } from '../../../modules/local/eggnogmapper/emapper'

include { INTERPROSCAN5_DOWNLOADDB                     } from '../../../modules/local/interproscan5/download_db'
include { INTERPROSCAN5_INTERPROSCAN as INTERPROSCAN5  } from '../../../modules/local/interproscan5/interproscan'

include { COMPLEMENT_GFF_WITH_INTERPROSCAN_GFF         } from '../../../modules/local/complement_gff_with_interproscan_gff'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_input
    functional_annotators
    interproscan5_db
    interproscan5_db_url

    main:

    if ( "eggnogmapper" in functional_annotators ) {

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // DOWNLOAD EGGNOG DB
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_eggnog_db = EGGNOGMAPPER_DOWNLOADDB( )

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // RUN EGGNOG MAPPER
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_eggnog_mapper_out = EGGNOGMAPPER_EMAPPER(
            ch_input.map { rec -> record(id: rec.id, fasta: rec.proteome, gff: rec.gff) },
            ch_eggnog_db
        )

        ch_input = ch_input.join( ch_eggnog_mapper_out, by: 'id' )

        // NOTE: the 'gff' key now holds the gff decorated by eggnog mapper
    }

    if ( "interproscan5" in functional_annotators ) {

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // IF INTERPROSCAN DB WAS NOT PROVIDED, DOWNLOADING IT
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        if ( interproscan5_db ) {
            interproscan_db = channel.fromPath( interproscan5_db, checkExists: true )
        } else {
            interproscan_db = INTERPROSCAN5_DOWNLOADDB( channel.value( interproscan5_db_url ) )
        }

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // RUNNING INTERPROSCAN 5
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_interproscan_out = INTERPROSCAN5( 
            ch_input, 
            interproscan_db
        )  

        ch_input = ch_input.join( 
            ch_interproscan_out.map { rec -> record(id: rec.id, interproscan_gff: rec.gff) }, 
            by: 'id' 
        )

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // COMPLEMENT EXISTING GFF WITH OUTPUT FROM INTERPROSCAN
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_complemented = COMPLEMENT_GFF_WITH_INTERPROSCAN_GFF( ch_input )
        ch_input = ch_input.join( ch_complemented, by: 'id' )
    }


    emit:
    ch_input
}

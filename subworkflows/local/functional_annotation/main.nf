nextflow.enable.types = true

include { EGGNOGMAPPER_DOWNLOADDB                      } from '../../../modules/local/eggnogmapper/download_db'
include { EGGNOGMAPPER_EMAPPER                         } from '../../../modules/local/eggnogmapper/emapper'

include { INTERPROSCAN6                                } from '../../../modules/local/interproscan6'

include { INTERPROSCAN5_DOWNLOADDB                     } from '../../../modules/local/interproscan5/download_db'
include { INTERPROSCAN5_INTERPROSCAN as INTERPROSCAN5  } from '../../../modules/local/interproscan5/interproscan'

include { COMPLEMENT_GFF_WITH_INTERPROSCAN_GFF         } from '../../../modules/local/complement_gff_with_interproscan_gff'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    proteome: Path
    gff: Path
}

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_input: Channel<Input>
    functional_annotators: List<String>
    eggnog_mapper_mode: String
    interproscan5_db: String
    interproscan5_db_url: String

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // EGGNOG MAPPER
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( "eggnogmapper" in functional_annotators ) {

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // DOWNLOAD EGGNOG DB
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_eggnog_db = EGGNOGMAPPER_DOWNLOADDB( eggnog_mapper_mode )
    
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // RUN EGGNOG MAPPER
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ch_eggnog_mapper_out = EGGNOGMAPPER_EMAPPER(
            ch_input.map { rec -> record(id: rec.id, fasta: rec.proteome, gff: rec.gff, eggnog_mapper_db: rec.eggnog_mapper_db) },
            ch_eggnog_db,
            eggnog_mapper_mode
        )

        ch_input = ch_input.join( ch_eggnog_mapper_out, by: 'id' )

        // NOTE: the 'gff' key now holds the gff decorated by eggnog mapper
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // INTERPROSCAN
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if ( "interproscan6" in functional_annotators || "interproscan5" in functional_annotators ) {
    
        if ( "interproscan6" in functional_annotators ) {
    
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // RUNNING INTERPROSCAN 6
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            ch_interproscan_out = INTERPROSCAN6( 
                ch_input.map { rec -> record(id: rec.id, fasta: rec.proteome) }
            )  
        
    
        } else { // interproscan5
    
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // IF INTERPROSCAN DB WAS NOT PROVIDED, DOWNLOADING IT
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            if ( interproscan5_db ) {
                interproscan_db = channel.fromPath( interproscan5_db, checkExists: true )
            } else {
                interproscan_db = INTERPROSCAN5_DOWNLOADDB( channel.value( interproscan5_db_url ) )
            }
    
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // RUN INTERPROSCAN 5
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
            ch_interproscan_out = INTERPROSCAN5( 
                ch_input.map { rec -> record(id: rec.id, fasta: rec.proteome) }, 
                interproscan_db
            )  

        }

        if ( functional_annotators.size() > 1 ) {

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // COMPLEMENT EXISTING GFF WITH OUTPUT FROM INTERPROSCAN
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
            ch_input = ch_input.join( 
                ch_interproscan_out.map { rec -> record(id: rec.id, interproscan_gff: rec.gff) }, 
                by: 'id' 
            )
            
            ch_complemented = COMPLEMENT_GFF_WITH_INTERPROSCAN_GFF( ch_input )
            ch_input = ch_input.join( ch_complemented, by: 'id' )

        }
    }


    emit:
    ch_input
}

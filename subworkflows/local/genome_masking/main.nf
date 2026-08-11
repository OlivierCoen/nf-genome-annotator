nextflow.enable.types = true

include { RED_RED as RED                                            } from '../../../modules/local/red/red'

include { REPEATMODELER_BUILDDATABASE as BUILDDATABASE              } from '../../../modules/local/repeatmodeler/builddatabase'
include { REPEATMODELER_REPEATMODELER as REPEATMODELER              } from '../../../modules/local/repeatmodeler/repeatmodeler'
include { REPEATMASKER_REPEATMASKER   as REPEATMASKER               } from '../../../modules/local/repeatmasker/repeatmasker'

//include { EARLGREY_DOWNLOADDB                                       } from '../../../modules/local/earlgrey/download_db'
//include { EARLGREY_EARLGREY as EARLGREY                             } from '../../../modules/local/earlgrey/earlgrey'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Genome {
    id: String
    fasta: Path
}

workflow GENOME_MASKING {

    take:
    ch_input: Channel<Genome>
    genome_masker: String

    main:

    ch_cannot_be_masked = channel.empty()

    if ( genome_masker == "red" ) {

        ch_masked = RED( ch_input )

    } else if ( genome_masker == "repeatmasker" ) {

        ch_repeatmodeler_db = BUILDDATABASE( ch_input )

        ch_repeatmodeler_out = REPEATMODELER( ch_repeatmodeler_db )

        // SOMETIMES REPEAT MODELER DOES NOT FIND FAMILIES
        // THE GENOME SHOULD NOT BE MASKED IN SUCH CASES

        ch_input = ch_input.join( ch_repeatmodeler_out, by: 'id' )

        ch_can_be_masked    = ch_input.filter{ rec -> rec.lib != null }
        ch_cannot_be_masked = ch_input.filter{ rec -> rec.lib == null }

        ch_masked = REPEATMASKER( ch_can_be_masked )

    } else if ( genome_masker == "earlgrey" ) {

        //EARLGREY_DOWNLOADDB()

        //EARLGREY


    }

    ch_masked = ch_input
                .join( ch_masked, by: 'id' )
                .map{ rec ->
                        // setting the softmasked genome as the default fasta file from now on
                        rec = rec + record(unmasked_fasta: rec.fasta)
                        rec = rec + record(fasta: rec.softmasked)
                        // return the record without the 'softmasked' field
                        rec.subMap(rec.keySet() - ['softmasked'])
                    }

    // records in ch_masked may have additional fields (like repeats_gff)
    // but these fields should not be used inside the workflow
    emit:
    masked  = ch_masked.mix( ch_cannot_be_masked )

}

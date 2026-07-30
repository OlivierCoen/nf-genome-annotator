nextflow.enable.types = true

include { REPEATMODELER_BUILDDATABASE as BUILDDATABASE              } from '../../../modules/local/repeatmodeler/builddatabase'
include { REPEATMODELER_REPEATMODELER as REPEATMODELER              } from '../../../modules/local/repeatmodeler/repeatmodeler'
include { REPEATMASKER_REPEATMASKER as REPEATMASKER                 } from '../../../modules/local/repeatmasker/repeatmasker'

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

    main:

    BUILDDATABASE( ch_input )

    REPEATMODELER( BUILDDATABASE.out )

    // SOMETIMES REPEAT MODELER DOES NOT FIND FAMILIES
    // THE GENOME SHOULD NOT BE MASKED IN SUCH CASES

    ch_input = ch_input.join( REPEATMODELER.out, by: 'id' )

    ch_can_be_masked    = ch_input.filter{ rec -> rec.lib != null }
    ch_cannot_be_masked = ch_input.filter{ rec -> rec.lib == null }

    REPEATMASKER( ch_can_be_masked )

    ch_masked = ch_input
                   .join( REPEATMASKER.out, by: 'id' )
                   .map{ rec -> 
                        // setting the softmasked genome as the default fasta file from now on
                        rec.unmasked_fasta = rec.fasta
                        rec.fasta = rec.softmasked
                        // return the record without the 'softmasked' field
                        rec.subMap(rec.keySet() - ['softmasked'])
                    }
                    
    // records in ch_masked may have additional fields (like repeats_gff)
    // but these fields should not be used inside the workflow
    emit:
    masked  = ch_masked.mix( ch_cannot_be_masked ) 

}

include { REPEATMODELER_BUILDDATABASE as BUILDDATABASE              } from '../../../modules/local/repeatmodeler/builddatabase'
include { REPEATMODELER_REPEATMODELER as REPEATMODELER              } from '../../../modules/local/repeatmodeler/repeatmodeler'
include { REPEATMASKER_REPEATMASKER as REPEATMASKER                 } from '../../../modules/local/repeatmasker/repeatmasker'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow GENOME_MASKING {

    take:
    ch_genome

    main:

    BUILDDATABASE( ch_genome )

    REPEATMODELER( BUILDDATABASE.out.db )

    // SOMETIMES REPEAT MODELER DOES NOT FIND FAMILIES
    // SO THE GENOME SHOULD NOT BE MASKED

    ch_genome_lib = ch_genome
                        .join ( REPEATMODELER.out.fasta, remainder: true )
                        .branch {
                            meta, genome, lib ->
                                can_be_masked: lib != null
                                    [ meta, genome, lib ]
                                cannot_be_masked: lib == null
                                    [ meta, genome ]
                        }

    REPEATMASKER( ch_genome_lib.can_be_masked )

    ch_masked_genomes = ch_genome_lib.cannot_be_masked
                            .mix ( REPEATMASKER.out.masked )


    emit:
    masked_genome           = ch_masked_genomes

}

include { REPEATMODELER_BUILDDATABASE as BUILDDATABASE              } from '../../../modules/nf-core/repeatmodeler/builddatabase'
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

    ch_versions = Channel.empty()

    BUILDDATABASE( ch_genome )

    REPEATMODELER( BUILDDATABASE.out.db )

    // SOMETIMES REPEAT MODELER DOES NOT FIND FAMILIES
    // SO THE GENOME DOES NOT HAVE TO BE MASKED

    ch_genome
        .join ( REPEATMODELER.out.fasta, remainder: true )
        .branch {
            meta, genome, lib ->
                can_be_masked: lib != null
                    [ meta, genome, lib ]
                cannot_be_masked: lib == null
                    [ meta, genome ]
        }
        .set { ch_genome_lib }

    REPEATMASKER( ch_genome_lib.can_be_masked )

    ch_genome_lib.cannot_be_masked
        .mix ( REPEATMASKER.out.masked )
        .set { ch_masked_genomes }


    ch_versions
        .mix ( BUILDDATABASE.out.versions )
        .set { ch_versions }


    emit:
    masked_genome           = ch_masked_genomes
    versions                = ch_versions

}


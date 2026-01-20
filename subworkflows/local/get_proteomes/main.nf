include { AGAT_SPEXTRACTSEQUENCES                       } from '../../../modules/local/agat/spextractsequences'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow GET_PROTEOMES {

    take:
    ch_gffs
    ch_genome
    codon_usage_id

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // JOINING WITH GENOME WHILE MAINTAINING A DIFFERENCE FOR THE MAIN GFF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_gff_genome = ch_gffs.map { meta, file -> [ meta.id, meta, file ] }
                        .combine(
                            ch_genome.map { meta, genome -> [ meta.id, meta, genome ] },
                            by: 0
                        )
                        .map {
                            id, meta, gff, meta2, genome -> [ meta, gff, genome ]
                        }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GETTING PROTOME
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    AGAT_SPEXTRACTSEQUENCES (
        ch_gff_genome,
        codon_usage_id,
        []
    )

    emit:
    proteomes         = AGAT_SPEXTRACTSEQUENCES.out.proteins

}

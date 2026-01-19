include { AGAT_SPEXTRACTSEQUENCES                       } from '../../../modules/local/agat/spextractsequences'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow GET_PROTEOMES {

    take:
    ch_gff
    ch_intermediate_gffs
    ch_genome
    codon_usage_id

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // JOINING WITH GENOME WHILE MAINTAINING A DIFFERENCE FOR THE MAIN GFF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_gff_genome = ch_gff.map { meta, file -> [ meta.id, meta + [main_annotation: true], file ] }
                        .mix(
                            ch_intermediate_gffs.map { meta, file -> [ meta.id, meta + [main_annotation: false], file ] }
                        )
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

    ch_proteomes = AGAT_SPEXTRACTSEQUENCES.out.proteins

    ch_main_proteome = ch_proteomes
                            .filter {
                                meta, file -> meta.main_annotation == true
                            }

    emit:
    proteomes         = ch_proteomes
    main_proteome     = ch_main_proteome

}

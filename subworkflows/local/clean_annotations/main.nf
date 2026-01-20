include { AGAT_SPFIXFEATURESLOCATIONSDUPLICATED   as AGAT_FIX_DUPLICATIONS                             } from '../../../modules/local/agat/spfixfeatureslocationsduplicated'
include { AGAT_SPKEEPLONGESTISOFORM               as AGAT_KEEP_LONGEST_ISOFORM                         } from '../../../modules/local/agat/spkeeplongestisoform'
include { AGAT_SPFIXOVERLAPPINGGENES              as AGAT_FIX_OVERLAPPING_GENES                        } from '../../../modules/local/agat/spfixoverlappinggenes'
include { AGAT_SPFILTERINCOMPLETEGENECODINGMODELS as AGAT_FILTER_INCOMPLETE_GENE_CODING_MODELS         } from '../../../modules/local/agat/spfilterincompletegenecodingmodels'
include { AGAT_SPFIXCDSPHASES                     as AGAT_FIX_CDS_PHASES                               } from '../../../modules/local/agat/spfixcdsphases'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow CLEAN_ANNOTATIONS {

    take:
    ch_annotation
    ch_genome
    skip_gff_keep_longest_isoform
    skip_gff_fix_overlapping_genes
    skip_gff_filter_incomplete_gene_models
    skip_gff_fix_cds_phases

    main:

    ch_intermediate_gffs = channel.empty()

    // remove redundant entries and convert all GTFs / GFFs to GFFs
    AGAT_FIX_DUPLICATIONS ( ch_annotation, [] )
    ch_gff = AGAT_FIX_DUPLICATIONS.out.gff

    ch_intermediate_gffs = ch_intermediate_gffs
                                .mix( ch_annotation )
                                .mix( AGAT_FIX_DUPLICATIONS.out.gff )

    if ( !skip_gff_keep_longest_isoform ) {
        AGAT_KEEP_LONGEST_ISOFORM ( ch_gff, [] )
        ch_gff = AGAT_KEEP_LONGEST_ISOFORM.out.gff
        ch_intermediate_gffs = ch_intermediate_gffs.mix( AGAT_KEEP_LONGEST_ISOFORM.out.gff )
    }

    if ( !params.skip_gff_fix_overlapping_genes ) {
        AGAT_FIX_OVERLAPPING_GENES ( ch_gff, [] )
        ch_gff = AGAT_FIX_OVERLAPPING_GENES.out.gff
        ch_intermediate_gffs = ch_intermediate_gffs.mix( AGAT_FIX_OVERLAPPING_GENES.out.gff )
    }

    if ( !params.skip_gff_filter_incomplete_gene_models ) {
        AGAT_FILTER_INCOMPLETE_GENE_CODING_MODELS ( ch_gff.join( ch_genome ), [] )
        ch_gff = AGAT_FILTER_INCOMPLETE_GENE_CODING_MODELS.out.gff
        ch_intermediate_gffs = ch_intermediate_gffs.mix( AGAT_FILTER_INCOMPLETE_GENE_CODING_MODELS.out.gff )
    }

    if ( !params.skip_gff_fix_cds_phases ) {
        AGAT_FIX_CDS_PHASES ( ch_gff.join( ch_genome ), [] )
        ch_gff = AGAT_FIX_CDS_PHASES.out.gff
        ch_intermediate_gffs = ch_intermediate_gffs.mix( AGAT_FIX_CDS_PHASES.out.gff ) // keeping it in case other cleaning steps are added afterwards
    }

    // removing main GFF from intermediate GFFs
    ch_intermediate_gffs = ch_intermediate_gffs
                            .combine ( ch_gff, by: 0 )
                            .filter {
                                meta, intermediate_gff, main_gff -> intermediate_gff != main_gff
                            }
                            .map {
                                meta, intermediate_gff, main_gff -> [ meta, intermediate_gff ]
                            }

    emit:
    gff                     = ch_gff
    intermediate_gffs       = ch_intermediate_gffs

}

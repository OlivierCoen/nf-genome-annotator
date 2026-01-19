include { AGAT_SPFIXFEATURESLOCATIONSDUPLICATED   as FIX_DUPLICATIONS                             } from '../../../modules/local/agat/spfixfeatureslocationsduplicated'
include { AGAT_SPKEEPLONGESTISOFORM               as KEEP_LONGEST_ISOFORM                         } from '../../../modules/local/agat/spkeeplongestisoform'
include { AGAT_SPFIXOVERLAPPINGGENES              as FIX_OVERLAPPING_GENES                        } from '../../../modules/local/agat/spfixoverlappinggenes'
include { AGAT_SPFILTERINCOMPLETEGENECODINGMODELS as FILTER_INCOMPLETE_GENE_CODING_MODELS         } from '../../../modules/local/agat/spfilterincompletegenecodingmodels'
include { AGAT_SPFIXCDSPHASES                     as FIX_CDS_PHASES                               } from '../../../modules/local/agat/spfixcdsphases'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow CLEAN_ANNOTATIONS {

    take:
    ch_annotation
    ch_genome

    main:

    ch_intermediate_gffs = channel.empty()

    // remove redundant entries and convert all GTFs / GFFs to GFFs
    FIX_DUPLICATIONS ( ch_annotation, [] )
    ch_gff = FIX_DUPLICATIONS.out.gff

    if ( !params.skip_gff_keep_longest_isoform ) {
        KEEP_LONGEST_ISOFORM ( ch_gff, [] )
        ch_gff = KEEP_LONGEST_ISOFORM.out.gff
        ch_intermediate_gffs = ch_intermediate_gffs.mix( KEEP_LONGEST_ISOFORM.out.gff )
    }

    if ( !params.skip_gff_fix_overlapping_genes ) {
        FIX_OVERLAPPING_GENES ( ch_gff, [] )
        ch_gff = FIX_OVERLAPPING_GENES.out.gff
        ch_intermediate_gffs = ch_intermediate_gffs.mix( FIX_OVERLAPPING_GENES.out.gff )
    }

    if ( !params.skip_gff_filter_incomplete_gene_models ) {
        FILTER_INCOMPLETE_GENE_CODING_MODELS ( ch_gff.join( ch_genome ), [] )
        ch_gff = FILTER_INCOMPLETE_GENE_CODING_MODELS.out.gff
        ch_intermediate_gffs = ch_intermediate_gffs.mix( FILTER_INCOMPLETE_GENE_CODING_MODELS.out.gff )
    }

    if ( !params.skip_gff_fix_cds_phases ) {
        FIX_CDS_PHASES ( ch_gff.join( ch_genome ), [] )
        ch_gff = FIX_CDS_PHASES.out.gff
        ch_intermediate_gffs = ch_intermediate_gffs.mix( FIX_CDS_PHASES.out.gff )
    }

    // removing main GFF from intermediate GFFs
    ch_intermediate_gffs = ch_intermediate_gffs
                            .join ( ch_gff )
                            .filter { meta, intermediate_gff, main_gff -> intermediate_gff != main_gff }
                            .map { meta, intermediate_gff, main_gff -> [ meta, intermediate_gff ] }

    emit:
    gff                     = ch_gff
    intermediate_gffs       = ch_intermediate_gffs

}

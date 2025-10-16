include { AGAT_SQREMOVEREDUNDANTENTRIES           as REMOVE_REDUNDANT_ENTRIES                     } from '../../../modules/local/agat/sq_remove_redundant_entries'
include { AGAT_SPKEEPLONGESTISOFORM               as KEEP_LONGEST_ISOFORM                         } from '../../../modules/local/agat/spkeeplongestisoform'
include { AGAT_SPFIXOVERLAPPINGGENES              as FIX_OVERLAPPING_GENES                        } from '../../../modules/local/agat/spfixoverlappinggenes'
include { AGAT_SPFILTERINCOMPLETEGENECODINGMODELS as FILTER_INCOMPLETE_GENE_CODING_MODELS         } from '../../../modules/local/agat/spfilterincompletegenecodingmodels'
include { AGAT_SPFIXCDSPHASES                     as FIX_CDS_PHASES                               } from '../../../modules/local/agat/spfixcdsphases'
include { AGAT_SPMANAGEIDS                        as FIX_IDS                                      } from '../../../modules/local/agat/spmanageids'

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

    // remove redundant entries and convert all GTFs / GFFs to GFFs
    REMOVE_REDUNDANT_ENTRIES ( ch_annotation, [] )
    REMOVE_REDUNDANT_ENTRIES.out.gff.set { ch_gff }

    if ( !params.skip_gff_keep_longest_isoform ) {
        KEEP_LONGEST_ISOFORM ( ch_gff, [] )
        KEEP_LONGEST_ISOFORM.out.gff.set { ch_gff }
    }

    if ( !params.skip_gff_fix_overlapping_genes ) {
        FIX_OVERLAPPING_GENES ( ch_gff, [] )
        FIX_OVERLAPPING_GENES.out.gff.set { ch_gff }
    }

    if ( !params.skip_gff_filter_incomplete_gene_models ) {
        FILTER_INCOMPLETE_GENE_CODING_MODELS ( ch_gff.join( ch_genome ), [] )
        FILTER_INCOMPLETE_GENE_CODING_MODELS.out.gff.set { ch_gff }
    }

    if ( !params.skip_gff_fix_cds_phases ) {
        FIX_CDS_PHASES ( ch_gff.join( ch_genome ), [] )
        FIX_CDS_PHASES.out.gff.set { ch_gff }
    }

    if ( !params.skip_gff_fix_ids ) {
        FIX_IDS ( ch_gff, [] )
        FIX_IDS.out.gff.set { ch_gff }
    }


    emit:
    gff                     = ch_gff

}


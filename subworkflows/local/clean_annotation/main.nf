nextflow.enable.types = true

include { AGAT_CONVERTSPGXF2GXF                   as AGAT_CONVERT_TO_GFF                               } from '../../../modules/local/agat/convertspgxf2gxf'
include { AGAT_SPFIXFEATURESLOCATIONSDUPLICATED   as AGAT_FIX_FEATURE_LOCATIONS_DUPLICATIONS           } from '../../../modules/local/agat/spfixfeatureslocationsduplicated'
include { AGAT_SPFIXOVERLAPPINGGENES              as AGAT_FIX_OVERLAPPING_GENES                        } from '../../../modules/local/agat/spfixoverlappinggenes'
include { AGAT_SPFILTERINCOMPLETEGENECODINGMODELS as AGAT_FILTER_INCOMPLETE_GENE_CODING_MODELS         } from '../../../modules/local/agat/spfilterincompletegenecodingmodels'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    structural_annotation: Path
    fasta: Path
}

def addGFFToIntermediateAnnotations( ch_input ){
    return ch_input.map { rec -> 
        rec.intermediate_annotations = rec.intermediate_annotations + [rec.gff]
        rec
    }
}


workflow CLEAN_ANNOTATION {

    take:
    ch_input: Channel<Input>
    gff_fix_feature_locations_duplicated: Boolean
    skip_gff_fix_overlapping_genes: Boolean
    skip_gff_filter_incomplete_gene_models: Boolean

    main:

    // for each modification, the workflow stores the version of the anntoation that is going to be replaced
    // in a specific list of intermediate GFFs

    ch_input.map { rec -> rec + record(intermediate_annotations: [rec.structural_annotation]) }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MANDATORY CLEANUP
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // remove redundant entries and convert all GTFs / GFFs to GFFs
    ch_out = AGAT_CONVERT_TO_GFF ( 
        ch_input.map { rec -> record(id: rec.id, gxf: rec.structural_annotation)}
    )
    ch_input = ch_input.join( ch_out, by: 'id' ).view()

    // each record has now a gff entry, that will be used in the subsequent steps

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // OPTIONAL CLEANUPS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( gff_fix_feature_locations_duplicated ) {
        ch_input = addGFFToIntermediateAnnotations( ch_input )
        ch_out = AGAT_FIX_FEATURE_LOCATIONS_DUPLICATIONS( ch_input )
        ch_input = ch_input.join( ch_out, by: 'id' )
    }

    if ( !skip_gff_fix_overlapping_genes ) {
        ch_input = addGFFToIntermediateAnnotations( ch_input )
        ch_out = AGAT_FIX_OVERLAPPING_GENES( ch_input )
        ch_input = ch_input.join( ch_out, by: 'id' )
    }

    if ( !skip_gff_filter_incomplete_gene_models ) {
        ch_input = addGFFToIntermediateAnnotations( ch_input )
        ch_out = AGAT_FILTER_INCOMPLETE_GENE_CODING_MODELS( ch_input )
        ch_input = ch_input.join( ch_out, by: 'id' )
    }

    emit:
    cleaned = ch_input

}

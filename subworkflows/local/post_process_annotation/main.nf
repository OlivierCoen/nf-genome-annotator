nextflow.enable.types = true

include { COMPLEMENT_ANNOTATION                                         } from '../complement_annotation'
include { CLEAN_ANNOTATION                                              } from '../clean_annotation'
include { ALTERNATIVE_ANNOTATIONS                                       } from '../alternative_annotation'

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


workflow POST_PROCESS_ANNOTATION {

    take:
    ch_input: Channel<Input>
    complement_annotation: Boolean
    skip_gff_cleaning: Boolean
    skip_alternative_annotations: Boolean
    gff_fix_feature_locations_duplicated: Boolean
    gff_fix_overlapping_genes: Boolean
    gff_filter_incomplete_gene_models: Boolean

    main:

    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // COMPLEMENTATION OF ANNOTATION (WHEN NECESSARY)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // complementation can only be done using the new structural annotation

    if ( complement_annotation ) {
        ch_complemented = COMPLEMENT_ANNOTATION( ch_input )
        ch_input = ch_input.join( ch_complemented, by: 'id' )
    }

    // storing the provided gff (if any)
    // filtering to keep only records that have at least a structural annotation or a gff
    ch_input = ch_input
                .filter { rec -> rec.structural_annotation != null }
                .map { rec -> rec.gff ? rec + record(previous_annotation: rec.gff) : rec }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CLEANING OF GFF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !skip_gff_cleaning ) {
    
        ch_cleaned = CLEAN_ANNOTATION (
            ch_input,
            gff_fix_feature_locations_duplicated,
            gff_fix_overlapping_genes,
            gff_filter_incomplete_gene_models
        )
        ch_input = ch_input.join( ch_cleaned, by: 'id' )
    
        // NOTE: now the annotation is under the 'gff' key

    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAKE ALTERNATIVE ANNOTATIONS (LONGEST ISOFORMS ONLY, ...)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !skip_alternative_annotations ) {

        ch_alternative_annotations = ALTERNATIVE_ANNOTATIONS( ch_input )
        ch_input = ch_input.join( ch_alternative_annotations, by: 'id' )

    }


    emit:
    ch_input

}

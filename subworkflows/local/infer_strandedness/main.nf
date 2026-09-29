nextflow.enable.types = true

include { AGAT_CONVERTSPGFF2BED as CONVERT_GFF_TO_BED     } from '../../../modules/local/agat/convertspgffbed'
include { SEQTK_SAMPLE as SAMPLE_READS                    } from '../../../modules/local/seqtk/sample'
include { RSEQC_INFEREXPERIMENT                           } from '../../../modules/local/rseqc/inferexperiment'
include { INFER_STRAND_FROM_XS                            } from '../../../modules/local/strandedness/infer_from_xs'
include { HISAT2_ALIGN                                    } from '../../../modules/local/hisat2/align'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    read_id: String
    reads: List<Path>
    fasta: Path
    reference_gff: Path
    index: Path
}


workflow INFER_STRANDEDNESS {

    take:
    ch_input: Channel<Input>

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SUBSAMPLE READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def sample_size = 2e5 // samples 200000 reads
    ch_sampled = SAMPLE_READS( 
        ch_input,
        sample_size
    )
    ch_input = ch_input.join( ch_sampled, by: 'read_id' )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ALIGNED SAMPLED READS TO GENOME
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ch_aligned = HISAT2_ALIGN( ch_input )

    ch_input = ch_input.join( ch_aligned, by: 'read_id' )

    ch_input_has_gff = ch_input.filter { rec -> rec.reference_gff != null }
    ch_input_no_gff  = ch_input.filter { rec -> rec.reference_gff == null }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SAMPLES HAVING A REFERENCE GFF (MORE PRECISE INFERENCE OF STRANDEDNESS)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // CONVERT REFERENCE GFF (WHEN PROVIDED) TO BED

    ch_gff_converted = CONVERT_GFF_TO_BED( 
        ch_input_has_gff.map { rec -> record(id: rec.id, gff: rec.reference_gff) }
    )
    ch_input_with_bed = ch_input_has_gff.join(ch_gff_converted, by: 'id')

    // RUN RSEQC TO INFER MAPPING METADATA
    RSEQC_INFEREXPERIMENT( ch_input_with_bed )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SAMPLES WITHOUT A REFERENCE GFF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    INFER_STRAND_FROM_XS

    

    

    emit:
    ch_input

}

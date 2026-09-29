nextflow.enable.types = true

include { AGAT_CONVERTSPGFF2BED as CONVERT_GXF_TO_BED     } from '../../../modules/local/agat/convertspgff2bed'
include { SEQTK_SAMPLE as SAMPLE_READS                    } from '../../../modules/local/seqtk/sample'
include { RSEQC_INFEREXPERIMENT                           } from '../../../modules/local/rseqc/inferexperiment'
include { INFER_STRAND_FROM_XS                            } from '../../../modules/local/strandedness/infer_from_xs'
include { INFER_STRAND_FROM_RSEQC                         } from '../../../modules/local/strandedness/infer_from_rseqc'
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
    reference_gtf: Path?
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
        ch_input.map { rec -> record(id: rec.read_id, reads: rec.reads) },
        sample_size
    )
    ch_input = ch_input.join(
        ch_sampled.map { rec -> record(read_id: rec.id, reads: rec.reads) }, 
        by: 'read_id'
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ALIGNED SAMPLED READS TO GENOME
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ch_aligned = HISAT2_ALIGN( ch_input )

    ch_input = ch_input.join( ch_aligned, by: 'read_id' )

    ch_input_has_gtf = ch_input.filter { rec -> rec.reference_gtf != null }
    ch_input_no_gtf  = ch_input.filter { rec -> rec.reference_gtf == null }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SAMPLES HAVING A REFERENCE GFF (MORE PRECISE INFERENCE OF STRANDEDNESS)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // CONVERT REFERENCE GFF (WHEN PROVIDED) TO BED

    ch_gff_converted = CONVERT_GXF_TO_BED( 
        ch_input_has_gtf.map { rec -> record(id: rec.id, gxf: rec.reference_gtf) }.unique()
    )
    ch_input_with_bed = ch_input_has_gtf.join(ch_gff_converted, by: 'id')

    // RUN RSEQC TO INFER MAPPING METADATA
    ch_report = RSEQC_INFEREXPERIMENT( ch_input_with_bed )

    ch_rseqc_strandedness = INFER_STRAND_FROM_RSEQC( ch_report )

    ch_input_has_gtf = ch_input_has_gtf.join(ch_rseqc_strandedness, by: 'read_id')

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SAMPLES WITHOUT A REFERENCE GFF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ch_xs_strandedness = INFER_STRAND_FROM_XS( ch_input_no_gtf )

    ch_input_no_gtf = ch_input_no_gtf.join(ch_xs_strandedness, by: 'read_id')
    
    emit:
    ch_input_has_gtf.mix(ch_input_no_gtf)

}

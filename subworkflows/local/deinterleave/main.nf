nextflow.enable.types = true

include { DETECT_INTERLEAVING                      } from '../../../modules/local/detect_interleaving'
include { BBMAP_DEINTERLEAVE                        } from '../../../modules/local/bbmap/deinterleave'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Input {
    id: String
    read_id: String
    reads: List<Path>
}


workflow DEINTERLEAVE {

    take:
    ch_single_reads: Channel<Input>

    main:
    
    ch_interleaving_status = DETECT_INTERLEAVING( ch_single_reads )

    ch_single_reads = ch_single_reads.join(ch_interleaving_status, by: 'read_id')

    ch_interleaved_reads = ch_single_reads.filter { rec -> rec.interleaved == true }
    ch_single_end_reads  = ch_single_reads.filter { rec -> rec.interleaved == false }

    ch_bbmap_deinterleave_out = BBMAP_DEINTERLEAVE( ch_interleaved_reads )
    ch_deinterleaved = ch_single_end_reads.join(ch_bbmap_deinterleave_out, by: 'read_id')
    

    emit:
    ch_single_end_reads.mix( ch_deinterleaved )

}

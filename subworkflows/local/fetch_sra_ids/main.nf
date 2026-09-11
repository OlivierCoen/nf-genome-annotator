nextflow.enable.types = true

include { GET_SRA_METADATA                   } from '../../../modules/local/get_sra_metadata'
include { GET_RANDOM_SAMPLE                  } from '../../../modules/local/get_random_sample'

record Input {
    id: String
    taxid: String
}


workflow FETCH_SRA_IDS {

    take:
    ch_input: Channel<Input>
    nb_short_read_sra_datasets: Integer
    nb_long_read_sra_datasets: Integer
    sra_allow_single_end: Boolean
    sra_random_seed: Integer
    
    main:

    ch_sra_id_files = GET_SRA_METADATA( 
        ch_input.map { rec -> rec.taxid }.unique(),
        nb_short_read_sra_datasets,
        nb_long_read_sra_datasets,
        sra_allow_single_end
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SELECTING RANDOMLY A CERTAIN NUMBER OF RNASEQ DATASETS FROM NCBI SRA
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_short_read_files_to_sample = ch_sra_id_files
                                    .filter{ rec -> rec.short_read_sra_ids_file != null }
                                    .combine( channel.value( nb_short_read_sra_datasets ) )
                                    .map{ rec, nb_to_sample -> record(
                                        taxid: rec.taxid, 
                                        type: 'short_reads', 
                                        input_file: rec.short_read_sra_ids_file,
                                        nb_to_sample: nb_to_sample
                                    ) }
                                    
    ch_long_read_files_to_sample = ch_sra_id_files
                                    .filter{ rec -> rec.long_read_sra_ids_file != null }
                                    .combine( channel.value( nb_long_read_sra_datasets ) )
                                    .map{ rec, nb_to_sample -> record(
                                        taxid: rec.taxid, 
                                        type: 'long_reads', 
                                        input_file: rec.long_read_sra_ids_file,
                                        nb_to_sample: nb_to_sample
                                    ) }
                         
    ch_sampled_sra_id_files = GET_RANDOM_SAMPLE(
        ch_short_read_files_to_sample.mix( ch_long_read_files_to_sample ),
        sra_random_seed
    )

    // ------------------------------------------------------------------------------------
    // ARRANGING CHANNELS
    // ------------------------------------------------------------------------------------

    ch_short_read_sra_id_files = ch_sampled_sra_id_files.filter{ rec -> rec.type == 'short_reads' }
    ch_long_read_sra_id_files  = ch_sampled_sra_id_files.filter{ rec -> rec.type == 'long_reads'  }

    ch_short_read_sra_ids = ch_short_read_sra_id_files
                            .map { rec -> record(
                                taxid: rec.taxid, 
                                fetched_short_read_sra_ids: rec.sampled.splitText().collect{ s -> s.strip() } 
                            ) }

    ch_long_read_sra_ids = ch_long_read_sra_id_files
                            .map { rec -> record(
                                taxid: rec.taxid, 
                                fetched_long_read_sra_ids: rec.sampled.splitText().collect{ s -> s.strip() } 
                            ) }

    ch_input = ch_input
                .join( ch_short_read_sra_ids, by: 'taxid', remainder: true )
                .join( ch_long_read_sra_ids, by: 'taxid', remainder: true )

    emit:
    ch_input
}

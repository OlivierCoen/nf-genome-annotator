nextflow.enable.types = true

include { FETCH_SRA_IDS                                 } from '../fetch_sra_ids'
include { DOWNLOAD_READS                                } from '../download_reads'
include { SHORT_READ_PREPARATION                        } from '../short_read_preparation'
include { LONG_READ_PREPARATION                         } from '../long_read_preparation'
include { BAM_SORT_INDEX_STATS                          } from '../bam_sort_index_stats'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Read {
    id: String
    reads: List<Path>
}

record Input {
    id: String
    supplied_short_reads: List<Read>
    supplied_long_reads: List<Read>
    supplied_short_read_sra_ids: List<String>
    supplied_long_read_sra_ids: List<String>    
    fasta: Path
    gff: Path
}


workflow PREPARE_RNASEQ_DATA {

    take:
    ch_input: Channel<Input>
    skip_fastqc: Boolean
    skip_fastqc_raw: Boolean
    skip_fastqc_cleaned: Boolean
    skip_umi_extract: Boolean
    skip_short_read_cleaning: Boolean
    nb_short_read_sra_datasets: Integer
    nb_long_read_sra_datasets: Integer
    sra_max_size: String
    sra_allow_single_end: Boolean
    sra_random_seed: Integer
    short_read_mapper: String
    ignore_existing_gff_for_mapping: Boolean
    

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // WHEN NEEDED, FETCH SRA IDS CORRESPONDING TO THE PROVIDED SPECIES
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( params.fetch_sra_rnaseq ) {

        ch_sra_ids = FETCH_SRA_IDS( 
            ch_input,
            nb_short_read_sra_datasets,
            nb_long_read_sra_datasets,
            sra_max_size,
            sra_allow_single_end,
            sra_random_seed
        )
        ch_input = ch_input.join( ch_sra_ids, by: 'id', remainder: true )
        
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DOWNLOAD READS FROM SRA / ENA
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // merging supplied and fetched SRA IDs
    ch_input = ch_input
                .map { rec ->
                    rec + record(
                        short_read_sra_ids: rec.supplied_short_read_sra_ids + rec.fetched_short_read_sra_ids,
                        long_read_sra_ids: rec.supplied_long_read_sra_ids + rec.fetched_long_read_sra_ids
                    )
                }

    ch_downloaded_reads = DOWNLOAD_READS( ch_input ) 
    ch_input = ch_input.join( ch_downloaded_reads, by: 'id', remainder: true )
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // PREPARE RNASEQ DATA FOR STRUCTURAL ANNOTATION (CLEANING AND / OR MAPPING)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_input = ch_input.map{ rec ->
        def downloaded_short_reads = rec.downloaded_short_reads ?: []
        def downloaded_long_reads = rec.downloaded_long_reads ?: []
        rec + record(
            short_reads: rec.supplied_short_reads + downloaded_short_reads,
            long_reads: rec.supplied_long_reads + downloaded_long_reads
        )
    }

    ch_short_reads_mapped = SHORT_READ_PREPARATION(
        ch_input.filter{ rec -> rec.short_reads.size() > 0 }, // pass only samples for which there are short reads
        params.skip_fastqc,
        params.skip_fastqc_raw,
        params.skip_fastqc_cleaned,
        params.skip_umi_extract,
        params.skip_short_read_cleaning,
        params.short_read_mapper,
        params.ignore_existing_gff_for_mapping
    )
    ch_input = ch_input.join( ch_short_reads_mapped, by: 'id', remainder: true )

    
    ch_long_reads_mapped = LONG_READ_PREPARATION(
        ch_input.filter{ rec -> rec.long_reads.size() > 0 }, // pass only samples for which there are long reads
        params.skip_fastqc,
        params.skip_fastqc_raw,
        params.skip_fastqc_cleaned,
        params.skip_long_read_cleaning
    )

    ch_input = ch_input.join( ch_short_reads_mapped, by: 'id', remainder: true )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // SORT ALL BAMS (SUPPLIED + NEWLY PRODUCED) AND GET MAPPING STATS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_input = ch_input.map{ rec ->
        def new_rnaseq_bams = rec.new_rnaseq_bams ?: []
        rec + record(bams: rec.supplied_rnaseq_bams + new_rnaseq_bams)
    }
    
    ch_sorted_bam = BAM_SORT_INDEX_STATS(
        ch_input.filter { rec -> rec.bams.size() > 0 }
    )
    ch_input = ch_input.join( ch_sorted_bam, by: 'id', remainder: true )

    emit:
    ch_input

}

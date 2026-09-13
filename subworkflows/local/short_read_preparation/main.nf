nextflow.enable.types = true

//
// Read QC, UMI extraction and trimming
//
include { FASTQC as FASTQC_RAW                     } from '../../../modules/local/fastqc'
include { FASTQC as FASTQC_CLEANED                 } from '../../../modules/local/fastqc'
include { UMITOOLS_EXTRACT                         } from '../../../modules/local/umitools/extract'
include { FASTP                                    } from '../../../modules/local/fastp'
include { AGAT_CONVERTSPGFF2GTF  as CONVERT_TO_GTF } from '../../../modules/local/agat/convertspgff2gtf'

include { FASTQ_ALIGN_HISAT2                       } from '../fastq_align_hisat2'
include { FASTQ_ALIGN_STAR                         } from '../fastq_align_star'



workflow SHORT_READ_PREPARATION {

    take:
    ch_input
    skip_fastqc: Boolean
    skip_fastqc_raw: Boolean
    skip_fastqc_cleaned: Boolean
    skip_umi_extract: Boolean
    skip_short_read_cleaning: Boolean
    short_read_mapper: String
    ignore_existing_gff_for_mapping: Boolean

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CLEAN AND CONVERT GFF / GTF TO GTF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !ignore_existing_gff_for_mapping ) {

        ch_input_with_gff    = ch_input.filter{ rec -> rec.gff != null }
        ch_input_without_gff = ch_input.filter{ rec -> rec.gff == null }

        ch_converted = CONVERT_TO_GTF( ch_input_with_gff )

        ch_input_with_gff = ch_input_with_gff.join( ch_converted, by: 'id')
        ch_input = ch_input_without_gff.mix( ch_input_with_gff )

    } else {
        ch_input = ch_input.map{ rec -> rec + record(gtf: null) }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // EXTRACT READS 
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // make channel containing one pair of fastq files (or one single fastq file in case of single end)
    // per element
    ch_reads = ch_input
                .flatMap{ rec -> rec.reads_to_map.collect{ subrec ->
                    record(
                        sample_id: rec.id,
                        fasta: rec.fasta,
                        gtf: rec.gtf,
                        id: subrec.id,
                        reads: subrec.reads
                    ) }
                }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // QC ON RAW READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !skip_fastqc && !skip_fastqc_raw ) {
        FASTQC_RAW( ch_reads )
    }

    if ( !skip_umi_extract ) {
        UMITOOLS_EXTRACT( ch_reads )
        ch_reads = UMITOOLS_EXTRACT.out
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CLEANING (TRIMMING + ADAPTER REMOVAL) + QC ON CLEANED READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !skip_short_read_cleaning ) {

        ch_cleaned_reads = FASTP( ch_reads )
        ch_reads = ch_reads.join( ch_cleaned_reads, by: 'id' )

        if ( !skip_fastqc && !skip_fastqc_cleaned ) {
            FASTQC_CLEANED( ch_reads )
        }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (short_read_mapper == 'hisat2') {

        ch_mapped = FASTQ_ALIGN_HISAT2(
            ch_reads,
            ignore_existing_gff_for_mapping
        )

    } else if (short_read_mapper == 'star') {

        ch_mapped = FASTQ_ALIGN_STAR(
            ch_reads,
            ignore_existing_gff_for_mapping
        )

    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GROUP BAM FILES BY SAMPLE ID
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_mapped = ch_mapped
                .map { rec -> tuple(rec.sample_id, rec.bam) }
                .groupTuple()
                .map { id, bams -> record(id: id, new_rnaseq_bams: bams) }

    

    emit:
    ch_mapped
}

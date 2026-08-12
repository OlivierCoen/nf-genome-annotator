nextflow.enable.types = true

include { FASTQ_FASTQC_UMITOOLS_FASTP    } from '../fastq_fastqc_umitools_fastp'
include { FASTQ_ALIGN_HISAT2             } from '../fastq_align_hisat2'
include { FASTQ_ALIGN_STAR               } from '../fastq_align_star'

include { AGAT_CONVERTSPGFF2GTF  as CONVERT_TO_GTF } from '../../../modules/local/agat/convertspgff2gtf'


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
    all_reads: List<Read>
    fasta: Path
    gff: Path
}


workflow MAP_RNASEQ_READS {

    take:
    ch_input: Channel<Input>
    skip_fastqc
    skip_umi_extract
    skip_trimming
    rnaseq_mapper
    ignore_existing_gff_for_mapping

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CLEAN AND CONVERT GFF / GTF TO GTF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !ignore_existing_gff_for_mapping ) {

        ch_input_with_gff    = ch_input.filter{ rec -> rec.gff != null }
        ch_input_without_gff = ch_input.filter{ rec -> rec.gff == null }

        CONVERT_TO_GTF( ch_input_with_gff )

        ch_input_with_gff = ch_input_with_gff.join( CONVERT_TO_GTF.out, by: 'id')
        ch_input = ch_input_without_gff.mix( ch_input_with_gff )

    } else {
        ch_input = ch_input.map{ rec -> rec + record(gtf: null) }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FASTQC & FASTP
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // make channel containing one pair of fastq files (or one single fastq file in case of single end)
    // per element
    ch_reads_to_map = ch_input
                        .flatMap{ rec -> rec.reads_to_map.collect{ subrec ->
                            record(
                                sample_id: rec.id,
                                fasta: rec.fasta,
                                gtf: rec.gtf,
                                id: subrec.id,
                                reads: subrec.reads
                            ) }
                        }

    FASTQ_FASTQC_UMITOOLS_FASTP(
        ch_reads_to_map,
        skip_fastqc,
        skip_umi_extract,
        skip_trimming
    )

    ch_prepared_reads = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (rnaseq_mapper == 'hisat2') {

        FASTQ_ALIGN_HISAT2(
            ch_prepared_reads,
            ignore_existing_gff_for_mapping
        )
        ch_mapped = FASTQ_ALIGN_HISAT2.out.mapped

    } else if (rnaseq_mapper == 'star') {

        FASTQ_ALIGN_STAR(
            ch_prepared_reads,
            ignore_existing_gff_for_mapping
        )
        ch_mapped = FASTQ_ALIGN_STAR.out.mapped

    }

    ch_mapped = ch_mapped
                .map { rec -> tuple(rec.sample_id, rec.bam) }
                .groupTuple()
                .map { id, bams -> record(id: id, new_rnaseq_bams: bams) }

    emit:
    mapped = ch_mapped

}

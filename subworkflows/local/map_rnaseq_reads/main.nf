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

record MappingInput {
    id: String
    reads: List<Path>
    fasta: Path
    gff: Path
}


workflow MAP_RNASEQ_READS {

    take:
    ch_input: Channel<MappingInput>
    skip_fastqc
    skip_umi_extract
    skip_trimming
    rnaseq_mapper
    ignore_existing_gff_for_mapping

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FASTQC & FASTP
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    FASTQ_FASTQC_UMITOOLS_FASTP(
        ch_input.map{ rec -> rec.subMap(['id', 'reads']) },
        skip_fastqc,
        skip_umi_extract,
        skip_trimming
    )

    ch_input = ch_input.join( FASTQ_FASTQC_UMITOOLS_FASTP.out.reads, by: 'id')

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CLEAN AND CONVERT GFF / GTF TO GTF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !ignore_existing_gff_for_mapping ) {
        CONVERT_TO_GTF( ch_input )
        ch_input = ch_input.join( CONVERT_TO_GTF.out, by: 'id')
    } else {
        ch_input = ch_input.map{ rec -> rec + record(gtf: null) }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MAP READS
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (rnaseq_mapper == 'hisat2') {

        FASTQ_ALIGN_HISAT2(
            ch_input,
            ignore_existing_gff_for_mapping
        )
        ch_mapped = FASTQ_ALIGN_HISAT2.out.mapped

    } else if (rnaseq_mapper == 'star') {

        FASTQ_ALIGN_STAR(
            ch_input,
            ignore_existing_gff_for_mapping
        )
        ch_mapped = FASTQ_ALIGN_STAR.out.mapped

    }

    emit:
    mapped = ch_mapped

}

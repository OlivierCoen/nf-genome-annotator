include { HISAT2_EXTRACTSPLICESITES     } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { HISAT2_EXTRACTEXONS           } from '../../../modules/local/hisat2/extractexons'
include { HISAT2_BUILD                  } from '../../../modules/local/hisat2/build'
include { HISAT2_ALIGN                  } from '../../../modules/local/hisat2/align'


workflow FASTQ_ALIGN_HISAT2 {
    take:
    ch_genome
    ch_reads
    ch_gtf
    ignore_existing_gtf_for_mapping

    main:

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // INDEX GENOME FOR HISAT2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( !ignore_existing_gtf_for_mapping ) {

        ch_gtf_to_extract = ch_reads
                                .cross( ch_gtf ) { v -> v[0][0] } // match only on id
                                .map{ // [[meta, reads], [meta2, index]]
                                    read_part, gtf_part -> 
                                        def meta = gtf_part[0]
                                        def gtf = gtf_part[1]
                                        [ meta, gtf ]
                                }

        HISAT2_EXTRACTSPLICESITES( ch_gtf_to_extract )
        ch_splicessites = HISAT2_EXTRACTSPLICESITES.out.txt
        
        HISAT2_EXTRACTEXONS( ch_gtf_to_extract )
        ch_exons = HISAT2_EXTRACTEXONS.out.txt
        
    } else {
        ch_splicessites = channel.of( [] )
        ch_exons = channel.of( [] )
    }

    ch_hisat2_build_input = ch_genome
                                .join( ch_splicessites, remainder: true ) // gives objects even if no corresponding splice sites
                                .join( ch_exons, remainder: true ) // gives objects even if no corresponding exons
                                .filter{ meta, genome, splicesites, exons -> genome != null }
                                .map{
                                    meta, genome, splicesites, exons ->
                                        [ meta, genome, splicesites?: [], exons?: [] ]
                                }

    HISAT2_BUILD( ch_hisat2_build_input )
    ch_index = HISAT2_BUILD.out.index
    
    //
    // Map reads with HISAT2
    //

    ch_hisat2_input =  ch_reads
                        .cross( ch_index ) { v -> v[0][0] } // match only on id, ignore single_end
                        .map{ // [[meta, reads], [meta2, index]]
                            read_part, index_part -> 
                                def meta = read_part[0]
                                def reads = read_part[1]
                                def index = index_part[1]
                                [ meta, reads, index ]
                        }
    
    HISAT2_ALIGN( ch_hisat2_input )

    emit:
    bam                 = HISAT2_ALIGN.out.bam // channel: [ val(meta), path(bam) ]
}

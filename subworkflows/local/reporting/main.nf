include { MULTIQC                                                 } from '../../../modules/nf-core/multiqc'


include { paramsSummaryMap                                        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                    } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                  } from '../../nf-core/utils_nfcore_pipeline'
include { formatVersionsToYAML                                    } from '../utils_nfcore_genome_annotator_pipeline'
include { methodsDescriptionText                                  } from '../utils_nfcore_genome_annotator_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow REPORTING {

    take:
    ch_versions
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    // ------------------------------------------------------------------------------------
    // DATA
    // ------------------------------------------------------------------------------------

    ch_multiqc_files = channel.empty()
                        .mix( channel.topic('hisat2_summary') )
                        .mix( channel.topic('star_log_final') )
                        .mix( channel.topic('mqc_busco_short_summaries_txt') )
                        .mix( channel.topic('mqc_mrna_with_isoforms_gff_stats') )
                        .mix( channel.topic('mqc_rna_with_isoforms_gff_stats') )
                        .mix( channel.topic('mqc_transcript_with_isoforms_gff_stats') )
                        .mix( channel.topic('mqc_mrna_without_isoforms_gff_stats') )
                        .mix( channel.topic('mqc_rna_without_isoforms_gff_stats') )
                        .mix( channel.topic('mqc_transcript_without_isoforms_gff_stats') )

    // ------------------------------------------------------------------------------------
    // VERSIONS
    // ------------------------------------------------------------------------------------

    // Collate and save software versions obtained from topic channels
    // TODO: use the nf-core functions when they are adapted to channel topics

    // Collate and save software versions
    formatVersionsToYAML ( Channel.topic('versions') )
        .mix ( softwareVersionsToYAML( ch_versions ) ) // mix with versions obtained from emit outputs
        .collectFile(storeDir: "${outdir}/pipeline_info", name: 'software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }


    // ------------------------------------------------------------------------------------
    // CONFIG
    // ------------------------------------------------------------------------------------

    ch_multiqc_config = Channel.fromPath( "$projectDir/assets/multiqc_config.yml", checkIfExists: true )

    summary_params = paramsSummaryMap( workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value( paramsSummaryMultiqc(summary_params) )

    ch_multiqc_custom_config = multiqc_config ?
                                    Channel.fromPath(multiqc_config, checkIfExists: true) :
                                    Channel.empty()

    ch_multiqc_logo = multiqc_logo ?
                        Channel.fromPath(multiqc_logo, checkIfExists: true) :
                        Channel.empty()

    ch_multiqc_custom_methods_description = multiqc_methods_description ?
                                                file(multiqc_methods_description, checkIfExists: true) :
                                                file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    ch_methods_description = Channel.value( methodsDescriptionText(ch_multiqc_custom_methods_description) )

    // Adding metadata to MultiQC
    ch_multiqc_files = ch_multiqc_files
                            .mix( ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml') )
                            .mix( ch_collated_versions )
                            .mix( ch_methods_description.collectFile( name: 'methods_description_mqc.yaml', sort: true ) )
                            
    // ------------------------------------------------------------------------------------
    // ADDING KEY TO JOIN ON
    // ------------------------------------------------------------------------------------

    ch_multiqc_file_list = ch_multiqc_files
                            .mix( ch_collated_versions )
                            .mix(
                                ch_methods_description.collectFile(
                                    name: 'methods_description_mqc.yaml',
                                    sort: true
                                )
                            )
                            .flatten()
                            .toSortedList()
                            .map{ list -> [ [id: 'Final report'], list ] }

    ch_multiqc_config_list = ch_multiqc_config
                                .mix( ch_multiqc_custom_config )
                                .toSortedList()
                                .map{ list -> [ [id: 'Final report'], list ] }

    ch_multiqc_logo = ch_multiqc_logo.map{ file -> [ [id: 'Final report'], file ] }

    // ------------------------------------------------------------------------------------
    // MULTIQC
    // ------------------------------------------------------------------------------------

    ch_multiqc_input = ch_multiqc_file_list
                        .join( ch_multiqc_config_list )
                        .join( ch_multiqc_logo )
                        .map { meta, files, configs, logo -> [ meta, files, configs, logo , [], [] ] }
                        
    
    
    MULTIQC ( ch_multiqc_input )
    

    emit:
    //report = MULTIQC.out.report
    report = channel.empty()
}

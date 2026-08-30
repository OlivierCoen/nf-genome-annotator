nextflow.enable.types = true

include { MULTIQC                                } from '../../../modules/nf-core/multiqc'

include { methodsDescriptionText                 } from '../utils_nfcore_genomeannotator_pipeline'
include { paramsSummaryMultiqc                   } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                 } from '../../nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap                       } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow REPORTING {

    take:
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    // ------------------------------------------------------------------------------------
    // VERSIONS
    // ------------------------------------------------------------------------------------

    // all modules send their versions to the 'versions' topic channel
    ch_versions = channel.empty()
    
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
                            .mix(topic_versions_string)
                            .collectFile(
                                storeDir: "${outdir}/pipeline_info",
                                name: 'nf_core_'  +  'genomeannotator_software_'  + 'mqc_'  + 'versions.yml',
                                sort: true,
                                newLine: true
                            )

    // ------------------------------------------------------------------------------------
    // PREPARE MULTIQC INPUT
    // ------------------------------------------------------------------------------------

    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    ch_multiqc_custom_config = multiqc_config ?
        channel.fromPath(multiqc_config, checkIfExists: true) :
        channel.empty()

    ch_multiqc_config_list = ch_multiqc_config
                                .mix( ch_multiqc_custom_config )
                                .map { files -> files.toSorted() }

    ch_multiqc_logo          = multiqc_logo ?
        channel.fromPath(multiqc_logo, checkIfExists: true) :
        channel.of([])

    summary_params      = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))

    ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')

    ch_multiqc_custom_methods_description = multiqc_methods_description ?
        file(multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    ch_methods_description     = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    // ------------------------------------------------------------------------------------
    // DATA
    // ------------------------------------------------------------------------------------

    ch_multiqc_files = channel.topic('multiqc')

    ch_multiqc_file_list = ch_multiqc_files
                            .groupTuple()
                            .combine( ch_collated_versions )
                            .combine(
                                ch_methods_description.collectFile(
                                    name: 'methods_description_mqc.yaml',
                                    sort: true
                                )
                            )
                            .map { id, file_list -> 
                                println file_list
                                [ id, file_list.flatten().toSorted() ] } // flatten and sort

    // ------------------------------------------------------------------------------------
    // MULTIQC
    // ------------------------------------------------------------------------------------

    ch_multiqc_input = ch_multiqc_file_list
                        .combine( ch_multiqc_config_list )
                        .combine( ch_multiqc_logo )
                        .map { meta, files, configs, logo -> [ meta, files, configs, logo , [], [] ] }
                        
    MULTIQC ( ch_multiqc_input )
    
    emit:
    report = MULTIQC.out.report
}

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

    // Collate and save software versions
    //

    def topic_versions_string = channel.topic("versions")
                                .map { process, tool, version ->
                                    [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
                                }
                                .groupTuple()
                                .map { process, tool_versions ->
                                    tool_versions.unique().sort()
                                    "${process}:\n${tool_versions.join('\n')}"
                                }

    ch_collated_versions = topic_versions_string
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
                                .collect()
                                .map { file_list -> file_list.toSorted() }

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

    ch_multiqc_files = channel.empty()
                        .mix( channel.topic('fastqc_multiqc') )
                        .mix( channel.topic('fastp_multiqc') )
                        .mix( channel.topic('hisat2_multiqc') )
                        .mix( channel.topic('star_multiqc') )
                        .mix( channel.topic('samtools_stat_multiqc') )
                        .mix( channel.topic('samtools_idxstat_multiqc') )
                        .mix( channel.topic('samtools_flagstat_multiqc') )
                        .mix( channel.topic('agat_structural_annotation_stats_multiqc') )
                        .mix( channel.topic('agat_functional_annotation_stats_multiqc').flatMap{ id, files -> files.collect{ file -> [id, file] } } )
                        .mix( channel.topic('busco_multiqc') )
                        .view{ v -> "multiqc files $v"}

    ch_multiqc_file_list = ch_multiqc_files
                            .groupTuple()
                            .combine( ch_collated_versions )
                            .combine(
                                ch_methods_description.collectFile(
                                    name: 'methods_description_mqc.yaml',
                                    sort: true
                                )
                            )
                            .map { data, version_file, description_file -> 
                                def id = data[0]
                                def meta = [id: id]
                                def file_list = data[1]
                                def data_files = file_list + [version_file, description_file]
                                [ meta, data_files.flatten().toSorted() ] 
                            } // flatten and sort for reproducibility

    // ------------------------------------------------------------------------------------
    // MULTIQC
    // ------------------------------------------------------------------------------------

    ch_multiqc_input = ch_multiqc_file_list
                        .combine( ch_multiqc_config_list )
                        .combine( ch_multiqc_logo )
                        .map { data, configs, logo -> 
                            def meta = data[0]
                            def data_files = data[1]
                            [meta, data_files, configs, logo, [], []] 
                        }
                        
    MULTIQC ( ch_multiqc_input )
    
    emit:
    report = MULTIQC.out.report.map{ meta, report -> record(id: meta.id, multiqc_report: report) }
}

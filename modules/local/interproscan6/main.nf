nextflow.enable.types = true

include { removeKeyValueArgs; 
          buildCacheDirURI;
          getWorkflowConfigPaths;
          workflowHasResumeFlag } from '../../../subworkflows/local/utils_nfcore_genomeannotator_pipeline'

// hard-coded versions to ensure reproducibility
def getNextflowVersion(){
    return '26.10.5'
}

def getInterproscanVersion(){
    return '6.0.3.0'
}


// largely inspired from nf-cascade:
// https://github.com/mahesh-panchal/nf-cascade

process INTERPROSCAN6 {
    
    tag "$id"
    label 'process_high'

    input:
        record(
            id: String,     
            fasta: Path
        )

    output:
        record(
            id: id,
            gff: file('results/*.gff3', optional: true)
        )
    
    topic:
        tuple('interproscan', id, file('results/*.tsv', optional: true))  >> 'additional_results'
        tuple('interproscan', id, file('results/*.xml', optional: true))  >> 'additional_results'
        tuple('interproscan', id, file('results/*.json', optional: true)) >> 'additional_results'
        tuple("${task.process}", 'nextflow',     getNextflowVersion())  >> 'versions'
        tuple("${task.process}", 'interproscan', getInterproscanVersion()) >> 'versions'

    exec:
    def args    = task.ext.args ?: ''

    // ---------------------------------------------------
    // CACHE
    // ---------------------------------------------------
    
    // Set cache directory so workflow can `-resume`
    def cache_dir = buildCacheDirURI(task.process, [id: id, fasta: fasta])
    def cache_path = file(cache_dir)
    assert cache_path.mkdirs()

    // ---------------------------------------------------
    // ENV VARS
    // ---------------------------------------------------
    
    // NXF_* and TOWER_* env vars inherited from a Tower/Seqera Platform launch break the nested run
    // Excluded so the nested run falls back to its own defaults, except shared, namespaced
    // storage locations we still want it to reuse.
    def nxf_passthrough = [
        'NXF_HOME',
        'NXF_ASSETS',
        'NXF_TEMP',
        'NXF_PLUGINS_DIR',
        'NXF_SINGULARITY_LIBRARYDIR',
        'NXF_CONDA_CACHEDIR',
        'NXF_SINGULARITY_CACHEDIR',
        'NXF_CHARLIECLOUD_CACHEDIR',
        'NXF_SPACK_CACHEDIR',
    ]
    def child_env = System.getenv()
        .findAll { k, v -> !(k.startsWith('NXF') || k.startsWith('TOWER')) || k in nxf_passthrough }
        .collect { k, v -> "${k}=${v}" }

    // setting nextflow version in env vars
    child_env << "NXF_VER=${getNextflowVersion()}"

    // ---------------------------------------------------
    // PROFILES
    // ---------------------------------------------------

    def valid_profiles = workflow.profile.tokenize(',').findAll { !it.startsWith('test_') }

    // ---------------------------------------------------
    // CONFIGS
    // ---------------------------------------------------

    // searches if the parent nextflow run command had config files provided with '-c'
    def config_paths = getWorkflowConfigPaths('-c')
    def config_arg = config_paths ? "-c $config_paths": ""

    // searches if the parent nextflow run command had config files provided with '-C'
    def overriding_config_paths = getWorkflowConfigPaths('-C')
    def overriding_config_arg = overriding_config_paths ? "-C $overriding_config_paths": ""

    // ---------------------------------------------------
    // BUILD COMMAND
    // ---------------------------------------------------
    
    // drop --input/--datadir and their values from the user-provided arg list
    def ips_args = removeKeyValueArgs([args: args, excluded_args: ["--input", "--datadir"]])

    def resume_arg = workflowHasResumeFlag() ? '-resume' : ''

    // Construct nextflow command
    def nxf_cmd = [
        'nextflow',
            '-log .nextflow.log',
            'run',
            'ebi-pf-team/interproscan6',
            "-r ${getInterproscanVersion()}",
            "-profile ${valid_profiles.join(',')}",
            config_arg,
            overriding_config_arg,
            '-ansi-log false',
            resume_arg,
            "--input $fasta",
            "--datadir data",
            ips_args.join(' '),
            "--outdir ${task.workDir}/results",
    ].minus("").join(" ")
        
    // Copy command to shell script in work dir for reference/debugging.
    file("${task.workDir}/nf-cmd.sh").text = nxf_cmd

    // ---------------------------------------------------
    // RUN
    // ---------------------------------------------------
    
    // Run nextflow command locally in cache directory
    def process = nxf_cmd.execute(child_env, cache_path.toFile())
    
    // Print process output to stdout and stderr (uncommed the line below for debugging)
    //process.consumeProcessOutput(System.out, System.err)
    process.waitFor()

    // ---------------------------------------------------
    // CHECK RUN STATUS
    // ---------------------------------------------------
    
    // Copy nextflow log to work directory
    cache_path.resolve(".nextflow.log").copyTo("${task.workDir}/nextflow.log")
    assert process.exitValue() == 0: 
        """
        ============================================================
        PIPELINE FAILED: ${id}
        Exit Code: ${process.exitValue()}
        Pipeline Log: ${task.workDir}/nextflow.log
        ============================================================
        """.stripIndent()

    // Clean cache of failed tasks
    def clean_cmd = ["/usr/bin/env", "bash", "-c", "nextflow clean -f -before last && find work -type d -empty -delete"]
    def clean_process = clean_cmd.execute(child_env, cache_path.toFile())
    // clean_process.consumeProcessOutput(System.out, System.err)
    clean_process.waitFor()
    assert clean_process.exitValue() == 0: 
        """
        ============================================================
        CACHE CLEAN FAILED: ${id}
        Exit Code: ${clean_process.exitValue()}
        ============================================================
        """.stripIndent()
        
}

nextflow.enable.types = true

process EGGNOGMAPPER_DOWNLOADDB {
    label 'process_medium'

    //storeDir "${workflow.projectDir}/.nextflow/cache/eggnogmapper_db"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b514ff051f837a2f06ab62a9a16592b0a2171917f1c7b11bbfef84c632e99f3/data':
        'community.wave.seqera.io/library/aria2_httpx_pigz_tenacity:0cf226d5365c27ae' }"

    input:
        eggnog_mapper_mode: String

    output:
        file("data", type: 'dir')
        
    topic:
        tuple("${task.process}", 'aria2', eval("aria2c -v | head -1 | sed 's/aria2 version //g'"))     >> 'versions'
        tuple("${task.process}", 'pigz',  eval("pigz --version 2>&1 | sed 's/pigz //g'"))              >> 'versions'
        tuple("${task.process}", 'httpx', eval('python3 -c "import httpx; print(httpx.__version__)"')) >> 'versions'

    script:
    // hardcoding a storeDir mechanism
    def store_dir = file("${workflow.projectDir}/.nextflow/cache/eggnogmapper_db/")
    
    if ( store_dir.isDirectory() && store_dir.listDirectory().size() > 0 ) {
        """
        ln -s ${store_dir} data
        """
    } else {
        """
        mkdir data
    
        download_eggnog_data.v2.py \\
            --db ${eggnog_mapper_mode} \\
            --out data \\
            --ncpus $task.cpus

        #########################
        # storing downloaded data
        #########################

        mkdir -p ${store_dir}
        mv data/*  ${store_dir}
        rm -rf data
        ln -s ${store_dir} data
        """
    }
}

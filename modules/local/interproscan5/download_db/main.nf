nextflow.enable.types = true

process INTERPROSCAN5_DOWNLOADDB {

    label 'process_medium'
    tag "${db_url.tokenize('/')[-1] - '.tar.gz'}"

    //storeDir "${workflow.projectDir}/.nextflow/cache/interproscan5_db"

    errorStrategy {
        if (task.exitStatus == 100) {
            log.warn("md5 checksum failed for Interproscan DB URL ${db_url}.")
            return 'retry'
        }
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95c0d3d867f5bc805b926b08ee761a993b24062739743eb82cc56363e0f7817d/data':
        'community.wave.seqera.io/library/aria2:1.37.0--3a9ec328469995dd' }"

    input:
        db_url: String

    output:
        file("*/data", type: 'dir')

    topic:
        tuple("${task.process}", 'aria2', eval("aria2c -v | head -1 | sed 's/aria2 version //g'")) >> 'versions'

    script:
    def filename = db_url.tokenize("/")[-1]
    // hardcoding a storeDir mechanism
    def db_id = db_url.tokenize('/')[-1] - '.tar.gz'
    def store_dir = file("${workflow.projectDir}/.nextflow/cache/${db_id}/")
    if ( store_dir.isDirectory() && store_dir.listDirectory().size() > 0 ) {
        """
        ln -s ${store_dir} ${db_id}/data
        """
    } else {
        """
        aria2c \\
            -s ${task.cpus} \\
            -x ${task.cpus} \\
            --max-tries=10 \\
            --retry-wait=30 \\
            --timeout=60 \\
            "${db_url}"
    
        echo "Checking md5"
        aria2c -c "${db_url}.md5"
        md5sum -c --status ${filename}.md5 && echo "Checksum: OK" || exit 100
    
        echo "Extracting archive"
        tar -pxzf ${filename}
    
        echo "Deleting archive"
        rm ${filename} ${filename}.md5

        #########################
        # storing downloaded data
        #########################
        
        # searching for the data directory
        data_dir=""
        for folder in */; do
            if [ -d "\${folder}data" ]; then
                data_dir="\${folder}data"
                break
            fi
        done
        echo "data dir: \$data_dir"
        
        # if not found, raise an error
        if [ -z "\$data_dir" ]; then
            echo "Error: data directory not found"
            exit 1
        fi

        # moving directory and creating a symlink
        mkdir -p ${store_dir}
        mv \${data_dir}/* ${store_dir}
        rm -rf \$data_dir
        ln -s ${store_dir} \$data_dir
        """
    }


}

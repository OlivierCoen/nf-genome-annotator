nextflow.enable.types = true

process FAMDB_DOWNLOAD_DFAM {

    tag "$taxid"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d541c440f6e2784d52642407ef15b3b7fa16035473403a026caa9e4c8536fc1d/data':
        'community.wave.seqera.io/library/aria2_pigz_python_h5py:4019bf13e8c1b4ad' }"

    input:
        record(
            taxid: String,
            taxid_lineage: String
        )
        
    output:
        record(
            taxid: taxid,
            dfam_db: file('dfam', type: 'dir')
        )

    topic:
        tuple("${task.process}", 'python', eval("python3 --version | sed 's/Python //'") )            >> 'versions'
        tuple("${task.process}", 'polars', eval('python3 -c "import h5py; print(h5py.__version__)"')) >> 'versions'
        tuple("${task.process}", 'famdb',  eval("famdb.py --help | grep version | cut -d' ' -f5"))    >> 'versions'
        tuple("${task.process}", 'aria2',  eval("aria2c -v | head -1 | sed 's/aria2 version //g'"))   >> 'versions'
        tuple("${task.process}", 'pigz',   eval("pigz --version 2>&1 | sed 's/pigz //g'"))            >> 'versions'

    script:
    """
    mkdir dfam

    ##################################################
    # DOWNLOADING DFAM ROOT FILE
    ##################################################
    
    download_dfam_4.0.py \\
        --fetch-root \\
        --output-dir dfam \\
        --ncpus ${task.cpus}

    ##################################################
    # FETCHING DATA TO DOWNLOAD FOR THE SPECIES OF INTEREST
    ##################################################

    # looping through the whole lineage of taxids, starting from the species
    # as soon as we found something in the Dfam db for this node, we break
    
    IFS=',' read -ra taxids <<< "$taxid_lineage"
    
    for lineage_txid in "\${taxids[@]}"; do

        echo "Searching data in Dfam database for taxid \$lineage_txid"
        
        famdb.py \\
            -i dfam \\
            check "\$lineage_txid" \\
            > famdb_check.out 2>&1
    
        if grep -q "No species found for search term" famdb_check.out; then
            echo "No data found for taxid \$lineage_txid"
            rm famdb_check.out
            continue
        else
            echo "Found data for taxid \$lineage_txid"
            cat famdb_check.out
            break
        fi
    done

    ##################################################
    # DOWNLOADING DATA FOR THE SPECIES OF INTEREST
    ##################################################
    
    download_dfam_4.0.py \\
        --famdb-check-output famdb_check.out \\
        --output-dir dfam \\
        --ncpus ${task.cpus}
    """

}

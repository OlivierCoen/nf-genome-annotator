nextflow.enable.types = true

def getContainerOptions(taxid) {
    def mapping = "${workflow.projectDir}/.nextflow/cache/dfam/${taxid}:/opt/conda/share/famdb-3.0.0/Libraries/"
    if ( workflow.containerEngine in ['singularity', 'apptainer'] ) { 
        return "-B $mapping"
    } else {
        return "-v $mapping"
    }
}

process EARLGREY_EARLGREY {
    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    //container "docker.io/tobybaril/earlgrey:latest-nodfam" // TODO: add specific version when available

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2d/2dc6599afbe7feb0cf30bdaf768a492c565cfa1f141042ba222d38a38ea8b340/data':
        '' }"


    //containerOptions "${getContainerOptions(taxid)}"
    
    input:
        record(
            id: String,
            taxid: String,
            fasta: Path,
            dfam_db: Path
        )

    output:
        record(
            id: id,
            softmasked: file("*.softmasked.fa.gz")
        )
  
    topic:
        tuple("${task.process}", 'earlgrey', eval("earlGrey | grep version | sed 's/earlGrey version //g'")) >> 'versions'

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "$id"
    def is_compressed = fasta.getExtension() == "gz"    ? true : false
    def fasta_name    = is_compressed                   ? fasta.getBaseName() : fasta.name
    """
    ################################################
    # CONFIGURATION OF REPEATMASKER.....
    ################################################
    # because RepeatMasker is old-fashioned and overly complicated to configure, we need to copy the script and its associated scripts to a writable location...
    # then and only then we can configure it and eventually use it...
    
    # 1 - Copying
    mkdir -p repeatmasker/bin repeatmasker/share/
    BIN_DIR=\$(dirname \$(which RepeatMasker))
    SHARE_DIR=\$(dirname \$BIN_DIR)/share
    REPEATMASKER_SHARE_DIR=\${SHARE_DIR}/RepeatMasker
    FAMDB_SHARE_DIR=\$(find \$SHARE_DIR -maxdepth 1 -type d -name "famdb-*" 2>/dev/null | sort -V | tail -n 1)
    
    cp \${BIN_DIR}/RepeatMasker repeatmasker/bin
    cp -r \$REPEATMASKER_SHARE_DIR repeatmasker/share/
    cp -r \$FAMDB_SHARE_DIR repeatmasker/share/

    # 2 - Symlink to the Dfam db
    FAMDB_DIRNAME=\$(basename \$FAMDB_SHARE_DIR)
    mkdir -p repeatmasker/share/\${FAMDB_DIRNAME}/Libraries
    cp -P ${dfam_db} repeatmasker/share/\${FAMDB_DIRNAME}/Libraries/famdb

    # 3 - Configuring RepeatMasker
    cd repeatmasker/share/RepeatMasker
    echo 'Y' | perl ./configure \\
        -trf_prgm \${BIN_DIR}/trf \\
        -rmblast_dir \$BIN_DIR \\
        -hmmer_dir \$BIN_DIR \\
        -default_search_engine rmblast
    cd ..; cd ..; cd ..

    # 4 - Changing PATH
    export PATH=\$PWD/repeatmasker/bin:\$PATH

    ################################################
    # UNGZIP FASTA IF NEEDED
    ################################################

    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    ################################################
    # RUN EARLGREY
    ################################################
    
    run_earlGrey.sh \\
        -g $fasta_name \\
        -o results \\
        -s ${prefix} \\
        -t ${task.cpus} \\
        -d yes \\
        ${args} \\
        || exit_code=\$?

    SOFTMASKED_OUTDIR="results/${prefix}_EarlGrey/${prefix}_summaryFiles"
    SOFTMASKED_OUTFILE=${prefix}.softmasked.fa
    
    if [ \$exit_code -eq 0 ]; then
        cp \${SOFTMASKED_OUTDIR}/*.fa */ \$SOFTMASKED_OUTFILE
    elif [ \$exit_code -eq 100 ]; then
        echo "Copying unmasked fasta to softmasked output"
        cp $fasta_name \$SOFTMASKED_OUTFILE
    else
        echo "earlGrey failed with exit code \$exit_code"
        exit \$exit_code
    fi

    echo "Compressing softmasked output"
    gzip \$SOFTMASKED_OUTFILE
    """

}

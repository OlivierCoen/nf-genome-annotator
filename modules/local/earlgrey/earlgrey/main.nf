nextflow.enable.types = true

process EARLGREY_EARLGREY {
    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1e39da20f65a43bc87cf59e52c1946929e616dcb07085201aca2cc5fa67ab91a/data':
        'community.wave.seqera.io/library/earlgrey_findutils_gzip_h5py_pruned:0e1759983126c6d4' }"

    
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
        tuple('earlgrey', id, file("results/*_EarlGrey/*_Database/*-families.fa",  optional: true)) >> 'additional_results'
        tuple('earlgrey', id, file("results/*_EarlGrey/*_Database/*-families.stk", optional: true)) >> 'additional_results'
        tuple('earlgrey', id, file("results/*_EarlGrey/*_Database/*-rmod.log",     optional: true)) >> 'additional_results'
        tuple('earlgrey', id, files("results/*_EarlGrey/*_summaryFiles/*",         optional: true)) >> 'additional_results'
        tuple("${task.process}", 'earlgrey', eval("earlGrey | grep version | sed 's/earlGrey version //g'")) >> 'versions'

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "$id"
    def is_compressed = fasta.getExtension() == "gz"    ? true : false
    def fasta_name    = is_compressed                   ? fasta.getBaseName() : fasta.name
    """
    ################################################
    # CONFIGURATION OF REPEATMODELER / REPEATMASKER.....
    ################################################
    # because RepeatMasker / RepeatMasker are really old-fashioned and overly complicated to configure
    # we need to copy the scripts and their helper scripts to a writable location...
    # then and only then we can configure it and eventually use it...

    WORKDIR=\$PWD
    
    # 1 - Copying binaries and software helper scripts / configurations
    BIN_DIR=\$(realpath "\${WORKDIR}/bin")
    SHARE_DIR=\$(realpath "\${WORKDIR}/share")
    mkdir \$BIN_DIR \$SHARE_DIR
    
    CONDA_BIN_DIR=\$(dirname \$(which RepeatMasker))
    CONDA_SHARE_DIR=\$(dirname \$CONDA_BIN_DIR)/share
   
    REPEATMODELER_CONDA_SHARE_DIR=\${CONDA_SHARE_DIR}/RepeatModeler
    REPEATMASKER_CONDA_SHARE_DIR=\${CONDA_SHARE_DIR}/RepeatMasker
    FAMDB_CONDA_SHARE_DIR=\$(find \$CONDA_SHARE_DIR -maxdepth 1 -type d -name "famdb-*" 2>/dev/null | sort -V | tail -n 1)
    FAMDB_DIRNAME=\$(basename \$FAMDB_CONDA_SHARE_DIR)

    # copying RepeatModeler / RepeatMasker / FamDB share directories (containing in particular the binaries) to the share directory
    cp -r \$REPEATMODELER_CONDA_SHARE_DIR \${SHARE_DIR}/
    cp -r \$REPEATMASKER_CONDA_SHARE_DIR \${SHARE_DIR}/
    cp -r \$FAMDB_CONDA_SHARE_DIR \${SHARE_DIR}/

    # 2 - Create newlinks to the real binaries in the new bin directory
    ln -s \${SHARE_DIR}/RepeatModeler/RepeatModeler \${BIN_DIR}/RepeatModeler
    ln -s \${SHARE_DIR}/RepeatModeler/BuildDatabase \${BIN_DIR}/BuildDatabase
    ln -s \${SHARE_DIR}/RepeatMasker/RepeatMasker \${BIN_DIR}/RepeatMasker
    
    # 3 - Symlink to the Dfam db
    mkdir -p \${SHARE_DIR}/RepeatMasker/Libraries
    cp -P ${dfam_db} \${SHARE_DIR}/RepeatMasker/Libraries/famdb

    # 4 - Setting path to FamDB in famdb.conf (for RepeatClassifier)
    sed -i "s|# FAMDB_DATA_DIR = /path/to/famdb/directory|FAMDB_DATA_DIR=\${SHARE_DIR}/RepeatMasker/Libraries/famdb|g" \${SHARE_DIR}/\${FAMDB_DIRNAME}/famdb.conf

    # 5 - Configuring RepeatMasker
    cd \${SHARE_DIR}/RepeatMasker
    echo 'Y' | perl ./configure \\
        -famdb_dir \${SHARE_DIR}/\${FAMDB_DIRNAME} \\
        -trf_prgm \${CONDA_BIN_DIR}/trf \\
        -rmblast_dir \$CONDA_BIN_DIR \\
        -hmmer_dir \$CONDA_BIN_DIR \\
        -default_search_engine rmblast
    cd \$WORKDIR

    # 6 - Configuring RepeatModeler
    cd \${SHARE_DIR}/RepeatModeler
    perl ./configure \\
        -repeatmasker_dir \${SHARE_DIR}/RepeatMasker \\
        -famdb_dir \${SHARE_DIR}/\${FAMDB_DIRNAME} \\
        -no_prompt
    cd \$WORKDIR

    # 7 - Modifying PATHS
    export PATH=\${BIN_DIR}:\${SHARE_DIR}/RepeatModeler:\${SHARE_DIR}/RepeatMasker:\${SHARE_DIR}/\${FAMDB_DIRNAME}:\$PATH
    export PERL5LIB="\${SHARE_DIR}/RepeatModeler:\${SHARE_DIR}/RepeatMasker:\${PERL5LIB:-}"
    
    ################################################
    # UNGZIP FASTA IF NEEDED
    ################################################

    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    ################################################
    # RUN EARLGREY
    ################################################

    exit_code=0
    
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
        if find \${SOFTMASKED_OUTDIR} -maxdepth 1 -name "*.softmasked.fasta" -print -quit | grep -q .; then
            echo "Found softmasked output, moving to \${SOFTMASKED_OUTFILE}"
            mv \${SOFTMASKED_OUTDIR}/*.softmasked.fasta \$SOFTMASKED_OUTFILE
        else
            echo "No softmasked output found, copying unmasked fasta to softmasked output"
            cp $fasta_name \$SOFTMASKED_OUTFILE
        fi
    elif [ \$exit_code -eq 100 ]; then
        echo "Copying unmasked fasta to softmasked output"
        cp $fasta_name \$SOFTMASKED_OUTFILE
    else
        echo "earlGrey failed with exit code \$exit_code"
        exit \$exit_code
    fi

    echo "Compressing softmasked output"
    gzip \$SOFTMASKED_OUTFILE

    rm -rf \$BIN_DIR \$SHARE_DIR
    """

}

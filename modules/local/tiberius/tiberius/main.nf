nextflow.enable.types = true

process TIBERIUS_TIBERIUS {
    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? ''
        : ''}"

    input:
        record(
            id: String,
            fasta: Path,
            tiberius_lineage: str,
            proteins,
            
        )
        
    
    topic:
        
        tuple("${task.process}", 'tiberius', '2.0.7') >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$id"
    def is_compressed = fasta.getExtension() == "gz"    ? true : false
    def fasta_name    = is_compressed                   ? fasta.getBaseName() : fasta.name
    """
    TIBERIUS_ENTRYPOINT=\$(find \${CONDA_PREFIX}/lib -name tiberius -type d | head -1)/main.py
    TIBERIUS_ROOT=\$(dirname \$TIBERIUS_ENTRYPOINT)
    MODEL_CFG=\$(dirname \${TIBERIUS_ROOT})/model_cfg/

    ################################################
    # DOWNLOAD MODEL CONFIG
    ################################################

    mkdir \$MODEL_CFG
    MODEL_URL="https://raw.githubusercontent.com/Gaius-Augustus/Tiberius/v2.0.7/model_cfg/${tiberius_lineage}.yaml"
    wget \$MODEL_URL -O \${MODEL_CFG}/${tiberius_lineage}.yaml

    ################################################
    # UNGZIP FASTA IF NEEDED
    ################################################

    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    ################################################
    # RUN
    ################################################

    \$TIBERIUS_ENTRYPOINT \\
        --genome $fasta_name \\
        --model_cfg $tiberius_lineage \\
        --out ${prefix}.tiberius.gff3 \\
        --work_dir \$PWD \\
        
        --threads ${task.cpus} \\
        --proteins [PROTEINS ...]
                            Protein FASTA input(s).
    
      --rnaseq_single [RNASEQ_SINGLE ...]
                            RNA-Seq single-end FASTQ input(s).
      --rnaseq_paired
         Paired-end RNA-Seq FASTQ input(s). From the CLI pass either a single quoted glob (e.g. --rnaseq_paired "RNA/*_{1,2}.fastq.gz") or exactly two FASTQ paths for one library (--rnaseq_paired r1.fq r2.fq). To list multiple explicit pairs use a params.yaml with a list of [r1, r2] pairs.
      --isoseq [ISOSEQ ...]
                              Iso-Seq FASTQ input(s).
        
    
    """

   
}

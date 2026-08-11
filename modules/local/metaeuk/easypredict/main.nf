nextflow.enable.types = true

process METAEUK_EASYPREDICT {
    tag "$id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a5/a5c4b28881614f8cb338a35e41bc5bf1ab7a77a0d94bfb5ad99c141af45400ce/data':
        'community.wave.seqera.io/library/metaeuk:7.bba0d80--479859525590824a' }"

    input:
        record(
            id: String,
            fasta: Path,
            db: Path
        )

    output:
        record(
            id: id,
            gff: file("*.gff")
        )

    topic:
        tuple('metaeuk', id, file("*.tsv"))       >> 'additional_results'
        tuple('metaeuk', id, file("*.codon.fas")) >> 'additional_results'
        tuple('metaeuk', id, file("*.fas"))       >> 'additional_results'
        tuple("${task.process}", 'metaeuk', eval("metaeuk | grep 'Version' | sed 's/metaeuk Version: //'"))  >> 'versions'

    script:
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "$id"
    """
    if [ -d ${db} ]; then
        ## if supplying an mmseqs database as a directory, metaeuk requires the basename of the database
        DBBASE=`find ${db}/ -name "*.version" -exec sh -c 'file=\$(basename {}); echo \${file%%.*}' \\;`
        DB=`echo "${db}/\${DBBASE}"`
    else
        DB=${db}
    fi

    metaeuk easy-predict \\
        ${fasta} \\
        \${DB} \\
        ${prefix} \\
        tmp/ \\
        ${args} \\
        --threads ${task.cpus}
    """
}

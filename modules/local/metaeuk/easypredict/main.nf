process METAEUK_EASYPREDICT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaeuk:6.a5d39d9--pl5321hf1761c0_2':
        'biocontainers/metaeuk:6.a5d39d9--pl5321hf1761c0_2' }"

    input:
    tuple val(meta), path(fasta), path(database)

    output:
    tuple val(meta), path("${prefix}.fas")      , emit: faa
    tuple val(meta), path("${prefix}.codon.fas"), emit: codon
    tuple val(meta), path("*.tsv")              , emit: tsv
    tuple val(meta), path("*.gff")              , emit: gff
    tuple val("${task.process}"), val('metaeuk'), eval("metaeuk | grep 'Version' | sed 's/metaeuk Version: //'"), topic: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    if [ -d ${database} ]; then
        ## if supplying an mmseqs database as a directory, metaeuk requires the basename of the database
        DBBASE=`find ${database}/ -name "*.version" -exec sh -c 'file=\$(basename {}); echo \${file%%.*}' \\;`
        DB=`echo "${database}/\${DBBASE}"`
    else
        DB=${database}
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

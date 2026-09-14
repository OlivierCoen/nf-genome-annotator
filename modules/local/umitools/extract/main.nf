nextflow.enable.types = true

process UMITOOLS_EXTRACT {
    tag "${id} :: ${read_id}"
    label "process_single"
    label "process_long"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/32/32476f0107d72dbd2210a4e56b2873abde07300025cc11052680475509d2db81/data' :
            'community.wave.seqera.io/library/umi_tools_future_matplotlib_numpy_pruned:1ee668bafc8c9f81' }"

    input:
        record(
            id: String, 
            read_id: String,
            reads: List<Path>
        )

    output:
        record(id: id, reads: file("*.fastq.gz"))


    topic:
        tuple('umitools_extract', id, file('*.log'))                                                     >> 'logs'
        tuple("${task.process}", 'umitools', eval("umi_tools --version | sed -n '/version:/s/.*: //p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$read_id"
    def single_end = reads.size() == 1 ? true : false
    if (single_end) {
        """
        umi_tools \\
            extract \\
            -I $reads \\
            -S ${prefix}.umi_extract.fastq.gz \\
            $args \\
            > ${prefix}.umi_extract.log
        """
    }  else {
        """
        umi_tools \\
            extract \\
            -I ${reads[0]} \\
            --read2-in=${reads[1]} \\
            -S ${prefix}.umi_extract_1.fastq.gz \\
            --read2-out=${prefix}.umi_extract_2.fastq.gz \\
            $args \\
            > ${prefix}.umi_extract.log
        """
    }
}

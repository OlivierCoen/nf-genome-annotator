nextflow.enable.types = true

process BRAKER3 {
    tag "$id"
    label 'process_high'

    // Re. Conda from the BRAKER team:
    // Warning: installing GeneMark-ETP for BRAKER in conda environments has lead to multiple problems reported by users (Issues!).
    // We can not offer support for conda installations. Please use the singularity image instead.
    // Two different images for short read / long read RNA-seq evidence
    container "${ long_read_rnaseq_bam ? 
        'docker://teambraker/braker3:isoseq' : 
        'docker.io/teambraker/braker3:v3.0.7.5' }"

    input:
        record(
            id: String,
            species: String,
            fasta: Path,
            proteins_fasta: Path?,
            short_read_rnaseq_bam: Path?,
            long_read_rnaseq_bam: Path?
        )

    stage:
        stageAs fasta, fasta.getExtension() == "gz" ? 'genome.fa.gz' : 'genome.fa'
        stageAs proteins_fasta, proteins_fasta.getExtension() == "gz" ? 'proteins.fa.gz' : 'proteins.fa'
        stageAs short_read_rnaseq_bam, 'short_read.rnaseq.bam'
        stageAs long_read_rnaseq_bam, 'long_read.rnaseq.bam'

    output:
        record(
            id: id,
            braker_gtf: file("workdir/braker.gtf", optional: true),
            braker_hintsfile: file("workdir/hintsfile.gff", optional: true),
        )

    topic:
        tuple('braker3', id, file('workdir/braker.log'))                       >> 'logs'
        tuple('braker3', id, file("workdir/braker.codingseq", optional: true)) >> 'additional_results'
        tuple('braker3', id, file("workdir/braker.aa", optional: true))        >> 'additional_results'
        tuple('braker3', id, file("workdir/braker.gff3", optional: true))      >> 'additional_results'
        tuple("${task.process}", 'braker3',   eval("braker.pl --version 2>/dev/null | sed 's/braker.pl version //'"))          >> 'versions'
        tuple("${task.process}", 'augustus',  eval("augustus --version |& sed -n 's/AUGUSTUS (\\(.*\\)) is a gene .*/\\1/p'")) >> 'versions'
        tuple("${task.process}", 'genemark',  eval("gmetp.pl | sed -n 's/ETP version \\(.*\\)/\\1/p'"))                        >> 'versions'
        tuple("${task.process}", 'prothint',  eval("prothint.py --version | sed 's/prothint.py //1'"))                         >> 'versions'

    script:
    def args               = task.ext.args                   ?: ''
    def prefix             = task.ext.prefix                 ?: "$id"
    // The number of CPUs cannot exceed 48, otherwise BRAKER warns that it could create problems with GeneMark
    def nb_threads         = Math.min(48, task.cpus)
    // FASTA
    def is_compressed      = fasta.getExtension() == "gz"    ? true : false
    def fasta_name         = is_compressed                   ? fasta.getBaseName() : fasta.name
    // PROTEINS
    def prot_is_compressed = proteins_fasta && proteins_fasta.getExtension() == "gz" ? true : false
    def prot_fasta_name    = proteins_fasta ? ( prot_is_compressed ? proteins_fasta.getBaseName() : proteins_fasta.name ) : null
    def prot_arg           = proteins_fasta ? "--prot_seq=$prot_fasta_name": ""
    // BAM
    // checking that either short read bam or isoseq bam was provided (or none), but not both
    def bam_arg = ''
    if ( short_read_rnaseq_bam && long_read_rnaseq_bam ) {
        error('Both short read and isoseq BAM files cannot be provided together')
    } else if (short_read_rnaseq_bam) {
        bam_arg = "--bam=$short_read_rnaseq_bam"
    } else if (long_read_rnaseq_bam) {
        bam_arg = "--bam=$long_read_rnaseq_bam"
    }
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    if [ -f $proteins_fasta -a "${prot_is_compressed}" == "true" ]; then
        gzip -c -d ${proteins_fasta} > ${prot_fasta_name}
    fi

    cp -r \$AUGUSTUS_CONFIG_PATH \\
        augustus_config

    chmod -R a+w \\
        augustus_config

    # keep only IDs in genome fasta header (and remove description)
    perl -p -e 's/^(>\\S+).*\$/\$1/' \\
        $fasta_name \\
        > ${prefix}.genome.masked.fasta

    braker.pl \\
        --genome ${prefix}.genome.masked.fasta \\
        --species "${species.replaceAll(/\s+/, "_")}_${id}" \\
        --workingdir workdir \\
        --AUGUSTUS_CONFIG_PATH "\$(pwd)/augustus_config" \\
        --threads $nb_threads \\
        $bam_arg \\
        $prot_arg \\
        $args
    """
}

process SRATOOLS_PREFETCH {
    tag "${meta.taxid} :: ${sra_id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/079d860bc0b82e189e4015772a7c6d8fe20d9b740058f854dd1513ce33460950/data' :
        'community.wave.seqera.io/library/sra-tools:3.2.1--2063130dadd340c5' }"

    input:
    tuple val(meta), val(sra_id)
    path ncbi_settings

    output:
    tuple val(meta), path("${two_first_letters}*", type: 'dir'),                                                                   emit: sra
    tuple val("${task.process}"), val('sratools'), eval("prefetch --version 2>&1 | grep -Eo '[0-9.]+'"),         topic: versions
    tuple val("${task.process}"), val('curl'), eval("curl --version | head -n 1 | sed 's/^curl //; s/ .*\$//'"), topic: versions

    script:
    args = task.ext.args ?: ''
    two_first_letters = sra_id[0..1]
    """
    export NCBI_SETTINGS="$PWD/!{ncbi_settings}"

    prefetch $args $sra_id

    # sometimes the SRA ID actually downloaded is different from the original one
    # but the two first letters stay the same
    DOWNLOADED_SRA=\$(find . -name "${two_first_letters}*" -type d)
    # check file integrity using vdb-validate or (when archive contains no checksums) md5sum
    vdb-validate \$DOWNLOADED_SRA > vdb-validate_result.txt 2>&1 || exit 100
    if grep -q "checksums missing" vdb-validate_result.txt; then
        VALID_MD5SUMS=\$(curl --silent --fail --location --retry 3 --retry-delay 60 "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?filetype=run&acc=\${DOWNLOADED_SRA}")
        FILES=\$(find \$DOWNLOADED_SRA -name "*" -type f)
        LOCAL_MD5SUMS=\$(md5sum \$FILES | cut -f1 -d' ')
        if ! grep -q -F -f <(echo "\$LOCAL_MD5SUMS") <(echo "\$VALID_MD5SUMS"); then
            echo "MD5 sum check failed" 1>&2
            exit 101
        fi
    fi
    """

    stub:
    """
    mkdir ${sra_id}
    touch ${sra_id}/${sra_id}.sra
    """
}

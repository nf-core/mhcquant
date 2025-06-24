process EASYPQP_CONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/easypqp:0.1.51--pyhdfd78af_0' :
        'biocontainers/easypqp:0.1.51--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(pepxml), path(spectra)
    path unimod

    output:
    tuple val(meta), path("*.psmpkl") , emit: psmpkl
    tuple val(meta), path("*.peakpkl"), emit: peakpkl
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    export MPLCONFIGDIR=/tmp/matplotlib
    export XDG_CACHE_HOME=/tmp/fontconfig-cache
    mkdir -p \$MPLCONFIGDIR \$XDG_CACHE_HOME

    easypqp convert \\
        --pepxml $pepxml \\
        --spectra $spectra \\
        --unimod $unimod \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        easypqp: \$(easypqp --version 2>&1 | grep -oP '(?<=easypqp, version )\\d+\\.\\d+\\.\\d+')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export MPLCONFIGDIR=/tmp/matplotlib
    export XDG_CACHE_HOME=/tmp/fontconfig-cache
    mkdir -p \$MPLCONFIGDIR \$XDG_CACHE_HOME

    touch "${prefix}.psmpkl"
    touch "${prefix}.peakpkl"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        easypqp: \$(easypqp --version 2>&1 | grep -oP '(?<=easypqp, version )\\d+\\.\\d+\\.\\d+')
    END_VERSIONS
    """
}

process MS2RESCORE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ms2rescore:3.1.5--pyh7e72e81_0':
        'biocontainers/ms2rescore:3.1.5--pyh7e72e81_0' }"

    // userEmulation settings when docker is specified
    containerOptions = (workflow.containerEngine == 'docker') ? '-u $(id -u) -e "HOME=${HOME}" -v /etc/passwd:/etc/passwd:ro -v /etc/shadow:/etc/shadow:ro -v /etc/group:/etc/group:ro -v $HOME:$HOME' : ''

    input:
    tuple val(meta), path(idxml), path(mzml), path(fasta)

    output:
    tuple val(meta), path("*ms2rescore.idXML") , emit: idxml
    tuple val(meta), path("*feature_names.tsv"), emit: feature_names
    tuple val(meta), path("*.html" )           , optional:true, emit: html
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_ms2rescore"

    """
    ms2rescore_cli.py \\
        --psm_file $idxml \\
        --spectrum_path . \\
        --output_path ${prefix}.idXML \\
        --processes $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MS²Rescore: \$(echo \$(ms2rescore --version 2>&1) | grep -oP 'MS²Rescore \\(v\\K[^\\)]+' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_ms2rescore"

    """
    touch ${prefix}.idXML
    touch ${meta.id}_feature_names.tsv
    touch ${meta.id}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MS²Rescore: \$(echo \$(ms2rescore --version 2>&1) | grep -oP 'MS²Rescore \\(v\\K[^\\)]+' ))
    END_VERSIONS
    """
}

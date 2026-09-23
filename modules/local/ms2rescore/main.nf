process MS2RESCORE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ms2rescore:4.0.2--pyhdfd78af_0':
        'biocontainers/ms2rescore:4.0.2--pyhdfd78af_0' }"

    // userEmulation settings when docker is specified
    containerOptions ( workflow.containerEngine == 'docker' ? '-u $(id -u) -e "HOME=${HOME}" -v /etc/passwd:/etc/passwd:ro -v /etc/shadow:/etc/shadow:ro -v /etc/group:/etc/group:ro -v $HOME:$HOME' : '' )

    input:
    tuple val(meta), path(idxml), path(mzml), path(fasta)

    output:
    tuple val(meta), path("*ms2rescore.idXML") , emit: idxml
    tuple val(meta), path("*feature_names.tsv"), emit: feature_names
    tuple val(meta), path("*.html" )           , optional:true, emit: html
    tuple val("${task.process}"), val('MS2Rescore'), eval('ms2rescore --version 2>/dev/null | tail -n1'), topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_ms2rescore"

    """
    touch ${prefix}.idXML
    touch ${meta.id}_feature_names.tsv
    touch ${meta.id}.html
    """
}

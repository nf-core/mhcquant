process SDRF_PIPELINES_PARSE_SDRF {
    label 'process_single'
    tag "${sdrf.baseName}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sdrf-pipelines:0.1.2--pyhdfd78af_0' :
        'quay.io/biocontainers/sdrf-pipelines:0.1.2--pyhdfd78af_0' }"

    input:
    path sdrf

    output:
    path "samplesheet.tsv"    , emit: samplesheet
    path "search_presets.tsv" , emit: search_presets
    tuple val("${task.process}"), val('sdrf-pipelines'), eval("parse_sdrf --version | cut -d ' ' -f 2"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parse_sdrf convert-mhcquant \\
        -s ${sdrf} \\
        -os samplesheet.tsv \\
        -op search_presets.tsv
    """

    stub:
    """
    touch samplesheet.tsv
    touch search_presets.tsv
    """
}

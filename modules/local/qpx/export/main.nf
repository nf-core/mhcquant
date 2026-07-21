process QPX_EXPORT {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qpx:1.0.2--pyhdfd78af_1' :
        'biocontainers/qpx:1.0.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(consensusxmls), path(sdrf), val(accession)

    output:
    path 'qpx/**'                                                    , emit: qpx
    tuple val("${task.process}"), val('qpx'), eval("qpxc --version"), topic: versions

    script:
    """
    openms_to_qpx.py \\
        --consensusxml ${consensusxmls.join(' ')} \\
        --accession ${accession} \\
        --outdir core/

    qpxc convert openms \\
        --qpx-dir core/ \\
        --sdrf-file ${sdrf} \\
        --project-accession ${accession} \\
        --output-folder qpx/

    qpxc validate --dataset-path qpx/
    """

    stub:
    """
    mkdir -p qpx
    touch qpx/${accession}.psm.parquet
    touch qpx/${accession}.feature.parquet
    touch qpx/${accession}.sdrf.tsv
    touch qpx/${accession}-project.json
    """
}

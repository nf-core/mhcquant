process QPX_EXPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qpx:1.0.2--pyhdfd78af_1' :
        'biocontainers/qpx:1.0.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(consensusxmls), path(sdrf), val(accession)

    output:
    path 'qpx/**'                                                    , emit: qpx
    tuple val("${task.process}"), val('qpx'), eval("qpxc --version | cut -d' ' -f2"), topic: versions

    when:
    task.ext.when == null || task.ext.when

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
        --output-prefix ${accession} \\
        --output-folder qpx/

    qpxc validate --dataset-path qpx/
    """

    stub:
    """
    mkdir -p qpx
    touch qpx/${accession}.{psm,feature,sample,run,ontology,provenance,dataset}.parquet
    """
}

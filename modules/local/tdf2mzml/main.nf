process TDF2MZML {
    tag "$meta.id"

    container "docker.io/mfreitas/tdf2mzml"

    input:
    tuple val(meta), path(tdf)

    output:
    tuple val(meta), path("*.mzML"), emit: mzml
    tuple val("${task.process}"), val('tdf2mzml'), val('0.4'), emit: versions_tdf2mzml, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${tdf.simpleName}"

    """
    tdf2mzml.py -i $tdf -o ${prefix}.mzML
    """

    stub:
    def prefix = task.ext.prefix ?: "${tdf.simpleName}"

    """
    touch ${prefix}.mzML
    """
}

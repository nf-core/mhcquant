process TDF2MZML {
    tag "$meta.id"

    container "docker.io/mfreitas/tdf2mzml:0.6.1_noentry"

    input:
    tuple val(meta), path(tdf)

    output:
    tuple val(meta), path("*.mzML"), emit: mzml
    tuple val("${task.process}"), val('python'), eval("python3 --version | cut -d ' ' -f2"), topic: versions
    tuple val("${task.process}"), val('tdf2mzml'), eval("tdf2mzml --version | cut -d' ' -f2"), topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "TDF2MZML module does not support Conda: tdf2mzml bundles the Bruker SDK and is only distributed as a Docker image. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${tdf.simpleName}"

    """
    tdf2mzml -i $tdf -o ${prefix}.mzML
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "TDF2MZML module does not support Conda: tdf2mzml bundles the Bruker SDK and is only distributed as a Docker image. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${tdf.simpleName}"

    """
    touch ${prefix}.mzML
    """
}

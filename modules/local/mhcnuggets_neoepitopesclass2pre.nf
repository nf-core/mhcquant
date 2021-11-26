// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MHCNUGGETS_NEOEPITOPESCLASS2RE {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0"
    } else {
        container "quay.io/biocontainers/mhcnuggets:2.3.2--py_0"
    }

    input:
        tuple val(meta), path(neoepitopes)

    output:
        tuple val(meta), path("*${prefix}*"), emit: preprocessed
        path "versions.yml"                 , emit: versions

    script:
        def prefix = options.suffix ? "${meta}_${options.suffix}" : "${meta}_mhcnuggets_preprocessed"

        """
        preprocess_neoepitopes_mhcnuggets.py --neoepitopes $neoepitopes --output ${prefix}
        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
        END_VERSIONS
        """
}

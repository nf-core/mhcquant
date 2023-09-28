process MHCNUGGETS_NEOEPITOPESCLASS2PRE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::mhcnuggets=2.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0' :
        'biocontainers/mhcnuggets:2.3.2--py_0' }"

    input:
        tuple val(meta), path(neoepitopes)

    output:
        tuple val(meta), path("*.csv")     , emit: preprocessed
        path "versions.yml"                , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${meta.id}_mhcnuggets_preprocessed"

        """
        preprocess_neoepitopes_mhcnuggets.py \\
            --neoepitopes $neoepitopes \\
            --output ${prefix}.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
        END_VERSIONS
        """
}

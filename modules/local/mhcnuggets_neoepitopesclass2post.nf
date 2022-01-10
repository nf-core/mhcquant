process MHCNUGGETS_NEOEPITOPESCLASS2POST {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0' :
        'quay.io/biocontainers/mhcnuggets:2.3.2--py_0' }"

    input:
        tuple val(meta), path(neoepitopes), path(predicted)

    output:
        tuple val(meta), path("*.csv"), emit: csv
        path "versions.yml"           , emit: versions

    script:

        """
        postprocess_neoepitopes_mhcnuggets.py --input $predicted --neoepitopes $neoepitopes

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
        END_VERSIONS
        """
}

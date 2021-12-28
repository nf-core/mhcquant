process MHCNUGGETS_PREDICTNEOEPITOPESCLASS2 {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0' :
        'quay.io/biocontainers/mhcnuggets:2.3.2--py_0' }"

    input:
        tuple val(meta), path(neoepitopes), val(alleles)

    output:
        tuple val(meta), path("*${prefix}"), emit: csv
        path "versions.yml"                , emit: versions

    script:
        def prefix           = task.ext.suffix ? "${meta}_${task.ext.suffix}" : "${meta}_predicted_neoepitopes_class_2"

        """
        mhcnuggets_predict_peptides.py --peptides $neoepitopes \\
            --alleles '$alleles' \\
            --output ${prefix}

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
        END_VERSIONS
        """
}

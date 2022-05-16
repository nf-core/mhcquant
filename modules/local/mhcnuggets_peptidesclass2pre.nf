process MHCNUGGETS_PEPTIDESCLASS2PRE {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0' :
        'quay.io/biocontainers/mhcnuggets:2.3.2--py_0' }"

    input:
        tuple val(meta), path(mztab)

    output:
        tuple val(meta), path("*_peptides")       , emit: preprocessed
        tuple val(meta), path('peptide_to_geneID'), emit: geneID
        path "versions.yml"                       , emit: versions

    script:
        def prefix           = task.ext.prefix ?: "${meta.sample}_preprocessed_mhcnuggets_peptides"

        """
        preprocess_peptides_mhcnuggets.py --mztab $mztab \\
            --output ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
        END_VERSIONS
        """
}

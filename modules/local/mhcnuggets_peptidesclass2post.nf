process MHCNUGGETS_PEPTIDESCLASS2POST {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::mhcnuggets=2.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0' :
        'biocontainers/mhcnuggets:2.3.2--py_0' }"

    input:
        tuple val(meta), path(peptides), path(peptide_to_geneID)

    output:
        tuple val(meta), path('*.csv'), emit: csv
        path "versions.yml"           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${meta.id}_postprocessed"

        """
        postprocess_peptides_mhcnuggets.py --input $peptides \\
            --peptides_seq_ID $peptide_to_geneID \\
            --output ${prefix}.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
        END_VERSIONS
        """
}

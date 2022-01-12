process MHCNUGGETS_PREDICTNEOEPITOPESCLASS2 {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0' :
        'quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0' }"

    input:
        tuple val(meta), path(neoepitopes), val(alleles)

    output:
        tuple val(meta), path("*.csv"), emit: csv
        path "versions.yml"           , emit: versions

    script:
        def prefix           = task.ext.prefix ?: "${meta}_predicted_neoepitopes_class_2"

        """
        mhcnuggets_predict_peptides.py --peptides $neoepitopes \\
            --alleles '$alleles' \\
            --output ${prefix}.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
        END_VERSIONS
        """
}

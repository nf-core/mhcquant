process MHCNUGGETS_PREDICTPEPTIDESCLASS2 {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0' :
        'quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0' }"

    input:
        tuple val(meta), path(peptides), val(alleles)

    output:
        tuple val(meta), path("*_class_2"), emit: csv
        path "versions.yml"               , emit: versions

    script:
        def prefix           = task.ext.suffix ? "${meta.sample}_${task.ext.suffix}_class_2" : "${meta.sample}_predicted_peptides_class_2"

        """
        mhcnuggets_predict_peptides.py --peptides $peptides \\
            --alleles '$alleles' \\
            --output ${prefix}

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            mhcflurry: \$(echo \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//') )
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
            fred2: \$(echo \$(python -c "import pkg_resources; print('fred2' + pkg_resources.get_distribution('Fred2').version)" | sed 's/^fred2//; s/ .*\$//'))
        END_VERSIONS
        """
}

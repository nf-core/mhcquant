process MHCFLURRY_PREDICTNEOEPITOPESCLASS1 {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0' :
        'quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0' }"

    input:
        tuple val(meta), val(allotypes), path(neoepitopes)

    output:
        tuple val(meta), path("*.csv"), emit: csv
        path "versions.yml"           , emit: versions

    script:
        def prefix           = task.ext.suffix ?: "${neoepitopes}_${meta}_predicted_neoepitopes_class_1"

        """
        mhcflurry-downloads --quiet fetch models_class1
        mhcflurry_neoepitope_binding_prediction.py '$allotypes' ${prefix}.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mhcflurry: \$(echo \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//') )
        END_VERSIONS
        """
}

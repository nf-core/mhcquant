process MHCFLURRY_PREDICTPEPTIDESCLASS1 {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::fred2=2.0.7 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c3f301504f7fa2e7bf81c3783de19a9990ea3001:12b1b9f040fd92a80629d58f8a558dde4820eb15-0' :
        'biocontainers/mulled-v2-c3f301504f7fa2e7bf81c3783de19a9990ea3001:12b1b9f040fd92a80629d58f8a558dde4820eb15-0' }"

    input:
        tuple val(meta), path(mztab), val(alleles)

    output:
        tuple val(meta), path("*.csv"), emit: csv
        path "versions.yml"           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.suffix ?: "${meta.id}_predicted_peptides_class_1"

        """
        mhcflurry-downloads --quiet fetch models_class1
        mhcflurry_predict_mztab.py '$alleles' $mztab ${prefix}.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
            mhcflurry: \$(echo \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//') )
            fred2: \$(echo \$(python -c "import pkg_resources; print('fred2' + pkg_resources.get_distribution('Fred2').version)" | sed 's/^fred2//; s/ .*\$//'))
        END_VERSIONS
        """
}

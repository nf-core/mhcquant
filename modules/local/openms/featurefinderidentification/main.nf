process OPENMS_FEATUREFINDERIDENTIFICATION  {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.4.1--h81ffffe_1' :
        'biocontainers/openms:3.4.1--h81ffffe_1' }"


    input:
    tuple val(meta), path(mzml), path(id_int), path(id_ext)

    output:
    tuple val(meta), path("*.featureXML"), emit: featurexml
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | grep -E '^Version' | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//'"), emit: versions_featurefinderidentification, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix    = task.ext.prefix ?: "${meta.id}_${meta.sample}_${meta.condition}"
    def quant_fdr = params.quantification_fdr ? "-id $id_int -id_ext $id_ext -svm:min_prob ${params.quantification_min_prob}" : "-id $id_ext"
    def args      = quant_fdr
    args          = args + (task.ext.args ? " ${task.ext.args}" : '')

    """
    FeatureFinderIdentification \\
        -in $mzml \\
        -out ${prefix}.featureXML \\
        -threads $task.cpus \\
        $args
    """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}_${meta.sample}_${meta.condition}"

        """
        touch ${prefix}.featureXML
        """
}

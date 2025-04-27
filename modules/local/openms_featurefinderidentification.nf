process OPENMS_FEATUREFINDERIDENTIFICATION  {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::openms=3.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.3.0--h0656172_8' :
        'biocontainers/openms:3.3.0--h0656172_8' }"


    input:
    tuple val(meta), path(mzml), path(id_int), path(id_ext)

    output:
    tuple val(meta), path("*.featureXML"), emit: featurexml
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix    = task.ext.prefix ?: "${meta.id}_${meta.sample}_${meta.condition}"
    def args      = task.ext.args  ?: ''
    def quant_fdr = params.quantification_fdr ? "-id $id_int -id_ext $id_ext -svm:min_prob ${params.quantification_min_prob}" : "-id $id_ext"
    args = args + " $quant_fdr"

    """
    FeatureFinderIdentification -in $mzml \\
        -out ${prefix}.featureXML \\
        -threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}_${meta.sample}_${meta.condition}"

        """
        touch ${prefix}.featureXML

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}

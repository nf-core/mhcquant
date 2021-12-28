process OPENMS_FEATUREFINDERIDENTIFICATION  {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::openms=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.6.0--h4afb90d_0' :
        'quay.io/biocontainers/openms:2.6.0--h4afb90d_0' }"

    input:
        tuple val(meta), path(id_quant_int), path(mzml), path(id_quant)

    output:
        tuple val(meta), path("*.featureXML"), emit: featurexml
        path "versions.yml"                  , emit: versions

    script:
        def prefix           = task.ext.suffix ? "${meta.sample}_${task.ext.suffix}" : "${meta.sample}_${meta.id}"

        if (!params.quantification_fdr){
            arguments = "-id $id_quant"
        } else {
            arguments = "-id $id_quant_int -id_ext $id_quant -svm:min_prob ${params.quantification_min_prob}"
        }

        """
        FeatureFinderIdentification -in $mzml \\
            -out ${prefix}.featureXML \\
            -threads $task.cpus \\
            ${arguments}

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}

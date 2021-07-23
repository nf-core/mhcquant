// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_FEATURELINKERUNLABELEDKD {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'RT_prediction', publish_id:'RT_prediction') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(meta), path(features)

    output:
        tuple val(meta), path("*.consensusXML"), emit: consensusxml   
        path  "*.version.txt", emit: version
        
    script:
        def software = getSoftwareName(task.process)
        def prefix = options.suffix ? "${meta.id}_${options.suffix}" : "${meta.id}_all_features_merged"

        """
            FeatureLinkerUnlabeledKD -in ${features} \\
                -out '${prefix}.consensusXML' \\
                -threads ${task.cpus}

            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt 
        
        """
}
// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process LINK_EXTRACTED_FEATURES {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'Intermediate_Results') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(features)

    output:
        tuple val("$Sample"), file("${Sample}_all_features_merged.consensusXML"), emit: consensusxml   
        path  "*.version.txt", emit: version
 
    when:
        !params.skip_quantification

    script:
    """
        FeatureLinkerUnlabeledKD -in ${features} \\
            -out '${Sample}_all_features_merged.consensusXML' \\
            -threads ${task.cpus}

        FileInfo --help &> openms.version.txt
    """
}
// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process LINK_EXTRACTED_FEATURES {
    publishDir "${params.outdir}/Intermediate_Results/"

    input:
    tuple val(id), val(Sample), val(Condition), file(features)

    output:
    tuple val("$Sample"), file("${Sample}_all_features_merged.consensusXML")
    
    when:
    !params.skip_quantification

    script:
    """
        FeatureLinkerUnlabeledKD -in ${features} \\
            -out '${Sample}_all_features_merged.consensusXML' \\
            -threads ${task.cpus}
    """

}
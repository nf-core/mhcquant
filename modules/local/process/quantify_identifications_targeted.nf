// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process QUANTIFY_IDENTIFICATION_TARGETED  {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(Sample), val(id), val(Condition), file(id_file_quant_int), file(mzml_quant), val(all_ids), file(id_file_quant)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_${id}.featureXML")

    when:
        !params.skip_quantification

    script:
    if (!params.quantification_fdr){
        """
            FeatureFinderIdentification -in ${mzml_quant} \\
                -id ${id_file_quant} \\
                -out ${Sample}_${id}.featureXML \\
                -threads ${task.cpus}
        """
    } else {
        """
            FeatureFinderIdentification -in ${mzml_quant} \\
                -id ${id_file_quant_int} \\
                -id_ext ${id_file_quant} \\
                -svm:min_prob ${params.quantification_min_prob} \\
                -out ${Sample}_${id}.featureXML \\
                -threads ${task.cpus}
        """   
    }
}
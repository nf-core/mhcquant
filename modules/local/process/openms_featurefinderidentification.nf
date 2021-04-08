// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process OPENMS_FEATUREFINDERIDENTIFICATION  {
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
        tuple val(Sample), val(id), val(Condition), file(id_file_quant_int), file(mzml_quant), val(all_ids), file(id_file_quant)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_${id}.featureXML"), emit: featurexml   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        
        if (!params.quantification_fdr){
            arguments = "-id ${id_file_quant}"
        } else {
            arguments = "-id ${id_file_quant_int} -id_ext ${id_file_quant} -svm:min_prob ${params.quantification_min_prob}"
        }

        """
            FeatureFinderIdentification -in ${mzml_quant} -out ${Sample}_${id}.featureXML -threads ${task.cpus} $arguments
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}   
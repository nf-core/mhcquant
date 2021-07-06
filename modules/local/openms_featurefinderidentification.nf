// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

def options    = initOptions(params.options)

//TODO: combine in a subflow --> when needs to be removed
process OPENMS_FEATUREFINDERIDENTIFICATION  {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'Intermediate_Results') }

    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(Sample), val(id), val(Condition), path(id_file_quant_int), path(mzml_quant), val(all_ids), path(id_file_quant)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), path("${Sample}_${id}.featureXML"), emit: featurexml   
        path  "*.version.txt", emit: version

    when:
        !params.skip_quantification

    script:
        def software = getSoftwareName(task.process)
        def prefix = options.suffix ? "${Sample}_${options.suffix}" : "${Sample}_${id}"
        
        if (!params.quantification_fdr){
            arguments = "-id ${id_file_quant}"
        } else {
            arguments = "-id ${id_file_quant_int} -id_ext ${id_file_quant} -svm:min_prob ${params.quantification_min_prob}"
        }

        """
            FeatureFinderIdentification -in ${mzml_quant} \\
                -out ${prefix}.featureXML \\
                -threads ${task.cpus} \\
                $arguments
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}   
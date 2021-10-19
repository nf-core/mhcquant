// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_FEATUREFINDERIDENTIFICATION  {
    tag "$meta.id"
    label 'process_intensive'

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
        tuple val(meta), path(id_quant_int), path(mzml), path(id_quant)

    output:
        tuple val(meta), path("*.featureXML"), emit: featurexml
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix = options.suffix ? "${meta.sample}_${options.suffix}" : "${meta.sample}_${meta.id}"

        if (!params.quantification_fdr){
            arguments = "-id ${id_quant}"
        } else {
            arguments = "-id ${id_quant_int} -id_ext ${id_quant} -svm:min_prob ${params.quantification_min_prob}"
        }

        """
            FeatureFinderIdentification -in ${mzml} \\
                -out ${prefix}.featureXML \\
                -threads ${task.cpus} \\
                $arguments
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}

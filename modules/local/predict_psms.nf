// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PREDICT_PSMS {
    tag "$meta"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'Intermediate_Results') }

    conda (params.enable_conda ? "bioconda::mhcflurry=2.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcflurry:2.0.1--pyh864c0ab_0"
    } else {
        container "quay.io/biocontainers/mhcflurry:2.0.1--pyh864c0ab_0"
    }

    input:
        tuple val(meta), path(perc_mztab), path(psm_mztab), val(allotypes)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path "versions.yml", emit: versions

    script:
        def prefix = options.suffix ? "${meta.id}_${options.suffix}" : "${meta.id}_peptide_filter"

    """
        mhcflurry-downloads --quiet fetch models_class1
        mhcflurry_predict_mztab_for_filtering.py ${params.subset_affinity_threshold} '${allotypes}' ${perc_mztab} ${psm_mztab} ${prefix}.idXML

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            mhcflurry: \$(mhcflurry-predict --version | sed 's/^mhcflurry //; s/ .*\$//' )
        END_VERSIONS
    """

}

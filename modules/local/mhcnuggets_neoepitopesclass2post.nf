// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MHCNUGGETS_NEOEPITOPESCLASS2POST {
    tag "$meta"
    label 'process_low'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'class_2_bindings', publish_id:'class_2_bindings') }

    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0"
    } else {
        container "quay.io/biocontainers/mhcnuggets:2.3.2--py_0"
    }


    input:
        tuple val(meta), path(neoepitopes), path(predicted)

    output:
        tuple val(meta), path("*.csv")  , emit: csv
        path "versions.yml"             , emit: versions

    script:

        """
            postprocess_neoepitopes_mhcnuggets.py --input ${predicted} --neoepitopes ${neoepitopes}
            cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
                mhcnuggets: \$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
            END_VERSIONS
        """
}

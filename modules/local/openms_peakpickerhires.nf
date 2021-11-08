// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_PEAKPICKERHIRES {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::openms=2.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.6.0--h4afb90d_0"
    } else {
        container "quay.io/biocontainers/openms:2.6.0--h4afb90d_0"
    }

    input:
        tuple val(meta), path(mzml)

    output:
        tuple val(meta), path("*.mzML"), emit: mzml
        path "versions.yml", emit: versions

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${mzml.baseName}_${options.suffix}" : "${mzml.baseName}"

        """
            PeakPickerHiRes -in ${mzml} \\
                -out ${prefix}.mzML \\
                -algorithm:ms_levels ${params.pick_ms_levels}

            cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
                openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
            END_VERSIONS
        """
}

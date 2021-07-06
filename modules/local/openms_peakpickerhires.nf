// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

process OPENMS_PEAKPICKERHIRES {
    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), path(mzml_file)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), path("*.mzML"), emit: mzml   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)

        """
            PeakPickerHiRes -in ${mzml_file} \\
                -out ${mzml_file.baseName}.mzML \\
                -algorithm:ms_levels ${params.pick_ms_levels}
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}
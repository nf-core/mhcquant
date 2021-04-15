// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process OPENMS_PEAKPICKERHIRES {
    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(mzml_file)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("*.mzML"), emit: mzml   
        path  "*.version.txt", emit: version

    when:
        params.run_centroidisation

    script:
        def software = getSoftwareName(task.process)

        """
            PeakPickerHiRes -in ${mzml_file} \\
                -out ${mzml_file.baseName}.mzML \\
                -algorithm:ms_levels ${params.pick_ms_levels}
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}
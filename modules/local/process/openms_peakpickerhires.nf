// Import generic module functions
include {  saveFiles; } from './functions'

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
        tuple val(id), val(Sample), val(Condition), file(mzml_unpicked)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${mzml_unpicked.baseName}.mzML"), emit: mzml   
        path  "*.version.txt", emit: version

    when:
        params.run_centroidisation

    script:
    """
        PeakPickerHiRes -in ${mzml_unpicked} \\
            -out ${mzml_unpicked.baseName}.mzML \\
            -algorithm:ms_levels ${params.pick_ms_levels}

        FileInfo --help &> openms.version.txt
    """
}
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_THERMORAWFILEPARSER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::thermorawfileparser::1.2.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/thermorawfileparser:1.2.3--1"
    } else {
        container "quay.io/biocontainers/thermorawfileparser:1.2.3--1"
    }

    input:
        tuple val(meta), path(rawfile)

    output:
        tuple val(meta), path("*.mzML"), emit: mzml
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${rawfile.baseName}_${options.suffix}" : "${rawfile.baseName}"

        """
            ThermoRawFileParser.sh -i=${rawfile} \\
                -f=2 \\
                -b=${prefix}.mzML
            ThermoRawFileParser.sh --version &> ThermoRawFileParser.version.txt
        """
}

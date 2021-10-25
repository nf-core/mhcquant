// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_THERMORAWFILEPARSER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::thermorawfileparser::1.2.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/thermorawfileparser:1.3.4--ha8f3691_0"
    } else {
        container "quay.io/biocontainers/thermorawfileparser:1.3.4--ha8f3691_0"
    }

    input:
        tuple val(meta), path(rawfile)

    output:
        tuple val(meta), path("*.mzML"), emit: mzml
        path "versions.yml", emit: versions

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${rawfile.baseName}_${options.suffix}" : "${rawfile.baseName}"

        """
            ThermoRawFileParser.sh -i=${rawfile} \\
                -f=2 \\
                -b=${prefix}.mzML
            > ThermoRawFileParser.version.txt

            cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
                thermorawfileparser: \$(ThermoRawFileParser.sh --version)
            END_VERSIONS
        """
}

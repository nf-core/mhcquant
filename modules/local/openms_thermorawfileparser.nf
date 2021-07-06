// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process OPENMS_THERMORAWFILEPARSER {

    conda (params.enable_conda ? "bioconda::thermorawfileparser::1.2.3" : null)

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/thermorawfileparser:1.2.3--1"
    } else {
        container "quay.io/biocontainers/thermorawfileparser:1.2.3--1"
    }

    input:
        tuple val(id), val(Sample), val(Condition), path(rawfile)
    
    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), path("*.mzML"), emit: mzml   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)

        """
            ThermoRawFileParser.sh -i=${rawfile} \\
                -f=2 \\
                -b=${rawfile.baseName}.mzML
                ThermoRawFileParser.sh &> ThermoRawFileParser.version.txt
        """

}
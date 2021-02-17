// Import generic module functions
include {  saveFiles; getSoftwareName } from './functions'

params.options = [:]

process RAW_FILE_CONVERSION {
    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(rawfile)
    
    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${rawfile.baseName}.mzML"), emit: mzml   
        path  "*.version.txt", emit: version

    script:
        """
        ThermoRawFileParser.sh -i=${rawfile} -f=2 -b=${rawfile.baseName}.mzML

        FileInfo --help &> openms.version.txt
        """
}
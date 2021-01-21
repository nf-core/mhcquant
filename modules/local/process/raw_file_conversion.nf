// Import generic module functions
include {  saveFiles } from './functions'

params.options = [:]

process RAW_FILE_CONVERSION {
    input:
    tuple val(id), val(Sample), val(Condition), file(rawfile)
    
    output:
    tuple val("$id"), val("$Sample"), val("$Condition"), file("${rawfile.baseName}.mzML")

    script:
    """
    ThermoRawFileParser.sh -i=${rawfile} -f=2 -b=${rawfile.baseName}.mzML
    """
}
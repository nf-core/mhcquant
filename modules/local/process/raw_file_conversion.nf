// Import generic module functions
include {  saveFiles } from './functions'

params.options = [:]

// TODO: include the right container... currently unknown
process RAW_FILE_CONVERSION {
    input:
    tuple val(id), val(Sample), val(Condition), file(rawfile)
    
    output:
    tuple val("$id"), val("$Sample"), val("$Condition"), file("${rawfile.baseName}.mzML")

    script:
    // TODO:(?) Change this to the FilConverted function of openMS
    """
    ThermoRawFileParser.sh -i=${rawfile} -f=2 -b=${rawfile.baseName}.mzML
    """
}
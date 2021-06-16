// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/*
 * Reformat design file and check validity
 */
process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path samplesheet
    
    output:
    path '*.csv'


    script:  
    
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv
    """
}

// Function to get list of [ meta, filenames ]
def get_samplesheet_paths(LinkedHashMap row) {
    def meta = [:]
    meta.id             = row.ID    

    def array = []
    if (!file(row.Filename).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> MS file does not exist!\n${row.Filename}"
    } else {
        array = [ meta, file(row.Filename) ]
    }
    return array    
}

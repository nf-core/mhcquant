//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:"\t" )
        .map { get_samplesheet_paths(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, filenames ]
def get_samplesheet_paths(LinkedHashMap row) {
    def meta = [:]
    meta.id             = row.ID
    meta.sample         = row.Sample
    meta.condition      = row.Condition
    meta.ext            = row.FileExt

    def array = []
    if (!file(row.Filename).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> MS file does not exist!\n${row.Filename}"
    } else {
        array = [ meta, file(row.Filename) ]
    }
    return array
}

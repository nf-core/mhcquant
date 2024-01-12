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
        .splitCsv ( header:true, sep:'\t' )
        .map { create_ms_channel(it) }
        .set { ms_runs }

    emit:
    ms_runs                                      // channel: [ val(meta), [ runs ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, filenames ]
def create_ms_channel(LinkedHashMap row) {
    def meta = [:]

    meta.id        = row.ID
    meta.sample    = row.Sample
    meta.condition = row.Condition
    meta.ext       = row.Extension

    // add path(s) of the data file(s) to the meta map
    def ms_file = file(row.ReplicateFileName)
    if (!ms_file.exists()) {
        exit 1, "ERROR: Please check input samplesheet -> MS file does not exist!\n${row.ReplicateFileName}"
    }

    meta.spectra   = ms_file.baseName
    ms_meta = [ meta, [ ms_file ] ]

    return ms_meta
}

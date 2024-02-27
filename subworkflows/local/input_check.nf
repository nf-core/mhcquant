//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    // Process the samplesheet and get the initial channels
    def initial_channels = SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv(header:true, sep:'\t')
        .map { row -> [row.Sample + '_' + row.Condition, row] }

    // Group by sample_condition and count
    def count_channels = initial_channels
        .groupTuple()
        .map { key, rows -> [rows, rows.size()] }
        .flatMap { rows, count -> rows.collect { [it, count] } }

    ms_runs = count_channels.map { row, group_count -> create_ms_channel(row, group_count) }

    emit:
    ms_runs                                      // channel: [ val(meta), [ runs ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, filenames ]
def create_ms_channel(LinkedHashMap row, int group_count) {
    def meta = [:]

    meta.id        = row.ID
    meta.sample    = row.Sample
    meta.condition = row.Condition
    meta.ext       = row.Extension
    meta.group_count = group_count

    // add path(s) of the data file(s) to the meta map
    def ms_file = file(row.ReplicateFileName)
    if (!ms_file.exists()) {
        exit 1, "ERROR: Please check input samplesheet -> MS file does not exist!\n${row.ReplicateFileName}"
    }

    meta.spectra   = ms_file.baseName
    ms_meta = [ meta, [ ms_file ] ]

    return ms_meta
}

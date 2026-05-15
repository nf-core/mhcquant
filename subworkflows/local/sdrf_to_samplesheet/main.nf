/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SDRF_TO_SAMPLESHEET: Convert SDRF to mhcquant samplesheet + download files from PRIDE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PARSESDRF_CONVERT     } from '../../../modules/nf-core/parsesdrf/convert/main'
include { PRIDEPY_FETCHSDRF     } from '../../../modules/nf-core/pridepy/fetchsdrf/main'
include { PRIDEPY_DOWNLOAD_FILE } from '../../../modules/local/pridepy/download_file/main'

workflow SDRF_TO_SAMPLESHEET {

    take:
    sdrf       // path: local SDRF file, or null
    pride_id   // val: PRIDE accession, or null

    main:

    // Resolve PRIDE accession for file downloads
    def resolved_accession = pride_id ?: inferPrideAccession(sdrf)

    // If pride_id given but no local SDRF, fetch from PRIDE
    if (pride_id && !sdrf) {
        PRIDEPY_FETCHSDRF(channel.of([[id: pride_id], pride_id]))
        ch_sdrf = PRIDEPY_FETCHSDRF.out.sdrf
    } else {
        ch_sdrf = channel.fromPath(sdrf, checkIfExists: true).map { f -> [[id: resolved_accession], f] }
    }

    // Convert SDRF to mhcquant samplesheet + search presets
    PARSESDRF_CONVERT(
        ch_sdrf.map { meta, sdrf_file -> [meta, sdrf_file, [], ''] },
        'mhcquant'
    )

    // Parse samplesheet to get file names for downloading
    ch_samplesheet_rows = PARSESDRF_CONVERT.out.mhcquant
        .map { _meta, samplesheet, _presets -> samplesheet }
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [
                id: row.ID as int,
                sample: row.Sample.toString(),
                condition: row.Condition.toString(),
                search_preset: row.SearchPreset
            ]
            [meta, row.ReplicateFileName]
        }

    // Download each file from PRIDE
    ch_to_download = ch_samplesheet_rows
        .map { meta, filename -> [meta, filename, resolved_accession] }

    PRIDEPY_DOWNLOAD_FILE(ch_to_download)

    // Write a validated samplesheet with local file paths
    ch_samplesheet_file = PRIDEPY_DOWNLOAD_FILE.out.downloaded_file
        .map { meta, downloaded_file ->
            [meta.id, meta.sample, meta.condition, downloaded_file, meta.search_preset].join('\t')
        }
        .collectFile(name: 'sdrf_samplesheet.tsv', seed: ['ID', 'Sample', 'Condition', 'ReplicateFileName', 'SearchPreset'].join('\t'), newLine: true)

    emit:
    samplesheet    = ch_samplesheet_file                                                       // path: samplesheet.tsv with local file paths
    search_presets = PARSESDRF_CONVERT.out.mhcquant.map { _meta, _samplesheet, presets -> presets }  // path: search_presets.tsv
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def inferPrideAccession(sdrf_path) {
    def name = file(sdrf_path).name
    def matcher = (name =~ /PXD\d{6,}/)
    if (matcher.find()) {
        return matcher.group()
    }
    error """\
        Could not infer PRIDE accession from SDRF filename: ${name}
        Please provide input as a PRIDE accession (e.g., --input PXD009752)
        or use an SDRF file named with the PXD accession (e.g., PXD009752.sdrf.tsv)
        """.stripIndent()
}

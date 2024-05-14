/*
 * Prepares the raw or compressed data holding spectra information for the subsequent database search.
 */

include { THERMORAWFILEPARSER    } from '../../modules/nf-core/thermorawfileparser/main'
include { UNTAR                  } from '../../modules/local/untar/main'
include { UNZIP                  } from '../../modules/local/unzip/main'
include { TDF2MZML               } from '../../modules/local/tdf2mzml'
include { GUNZIP                 } from '../../modules/nf-core/gunzip/main'
include { OPENMS_PEAKPICKERHIRES } from '../../modules/nf-core/openms/peakpickerhires/main'

workflow PREPARE_SPECTRA {
    take:
        ch_samplesheet

    main:
        ch_versions = Channel.empty()

        ch_samplesheet
        .branch {
            meta, filename ->
                raw : meta.ext == 'raw'
                    return [ meta, filename ]
                mzml : meta.ext == 'mzml'
                    return [ meta.subMap('id', 'sample', 'condition', 'group_count', 'spectra'), filename ]
                mzml_gz : meta.ext == 'mzML.gz'
                    return [ meta.subMap('id', 'sample', 'condition', 'group_count', 'spectra'), filename ]
                d : meta.ext == 'd'
                    return [ meta.subMap('id', 'sample', 'condition', 'group_count', 'spectra'), filename ]
                d_tar : meta.ext == 'd.tar' | meta.ext == 'd.tar.gz'
                    return [ meta.subMap('id', 'sample', 'condition', 'group_count', 'spectra'), filename ]
                d_zip : meta.ext == 'd.zip'
                    return [ meta.subMap('id', 'sample', 'condition', 'group_count', 'spectra'), filename ]
                other : true }
        .set { branched_ms_files }

        // Raw file conversion
        THERMORAWFILEPARSER(branched_ms_files.raw)
        ch_versions = ch_versions.mix(THERMORAWFILEPARSER.out.versions)

        // Decompress timsTOF archive for data conversion
        UNTAR(branched_ms_files.d_tar)
        ch_versions = ch_versions.mix(UNTAR.out.versions)

        UNZIP(branched_ms_files.d_zip)
        ch_versions = ch_versions.mix(UNZIP.out.versions)

        ch_tdf_files = branched_ms_files.d
                            .mix(UNTAR.out.untar,
                                UNZIP.out.unzipped_archive)

        // timsTOF data conversion
        TDF2MZML(ch_tdf_files)
        ch_versions = ch_versions.mix(TDF2MZML.out.versions)

        // Gunzip mzML files
        GUNZIP(branched_ms_files.mzml_gz)
        // Initialize channel for ms files that do not need to be converted
        ch_ms_files = branched_ms_files.mzml
                        .mix(GUNZIP.out.gunzip,
                            THERMORAWFILEPARSER.out.spectra,
                            TDF2MZML.out.mzml)

        // Optional: Run Peak Picking as Preprocessing
        if (params.run_centroidisation) {
            OPENMS_PEAKPICKERHIRES(ch_ms_files)
            ch_versions = ch_versions.mix(OPENMS_PEAKPICKERHIRES.out.versions)
            ch_mzml_file = OPENMS_PEAKPICKERHIRES.out.mzml
        } else {
            ch_mzml_file = ch_ms_files
        }

    emit:
        mzml = ch_mzml_file
        versions = ch_versions
}

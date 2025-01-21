/*
 * Generates spectrum library for DIA-based searches
 */

//
// MODULE: Loaded from modules/local/
//

include { EASYPQP_CONVERT                         } from '../../modules/local/easypqp/convert'
include { EASYPQP_LIBRARY;
        EASYPQP_LIBRARY as EASYPQP_LIBRARY_GLOBAL } from '../../modules/local/easypqp/library'

//
// MODULE: Installed directly from nf-core/modules
//

workflow SPECLIB {

    take:
        fdrfiltered_comet_idxml
        mzml

    main:
        ch_versions = Channel.empty()

    // Load unimod tables (Future:)
    unimod = file("$projectDir/assets/250120_unimod_tables.xml", checkIfExists: true)

    // Convert psms and spectra to pickle files
    EASYPQP_CONVERT(fdrfiltered_comet_idxml.join(mzml), unimod)
    ch_versions = ch_versions.mix(EASYPQP_CONVERT.out.versions)

    EASYPQP_CONVERT.out.psmpkl
        .map { meta, psmpkl -> [groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), psmpkl] }
        .groupTuple()
        .set { ch_psmpkl }
    EASYPQP_CONVERT.out.peakpkl
        .map { meta, peakpkl -> [groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), peakpkl] }
        .groupTuple()
        .set { ch_peakpkl }

    // Generate spectrum library for each sample-condition pair
    EASYPQP_LIBRARY(ch_psmpkl.join(ch_peakpkl))
    ch_versions = ch_versions.mix(EASYPQP_LIBRARY.out.versions)

    // Generate spectrum library for all MSruns in the samplesheet
    if (params.global_fdr) {
        EASYPQP_CONVERT.out.psmpkl
            .map { meta, psmpkl -> [[id: "global"], psmpkl] }
            .groupTuple()
            .set { ch_global_psmpkl }
        EASYPQP_CONVERT.out.peakpkl
            .map { meta, peakpkl -> [[id: "global"], peakpkl] }
            .groupTuple()
            .set { ch_global_peakpkl }
        EASYPQP_LIBRARY_GLOBAL(ch_global_psmpkl.join(ch_global_peakpkl))
    }

    emit:
        versions = ch_versions
}

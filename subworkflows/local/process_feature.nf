/*
 * Perform the quantification of the samples when the parameter --skip_quantification is not provided
 */

include { OPENMS_IDMERGER }                                                 from '../../modules/local/openms_idmerger'
include { OPENMS_FEATUREFINDERIDENTIFICATION }                              from '../../modules/local/openms_featurefinderidentification'
include { OPENMS_FEATURELINKERUNLABELEDKD }                                 from '../../modules/local/openms_featurelinkerunlabeledkd'
include { OPENMS_IDCONFLICTRESOLVER }                                       from '../../modules/local/openms_idconflictresolver'
include { OPENMS_TEXTEXPORTER as OPENMS_TEXTEXPORTER_QUANT }                from '../../modules/local/openms_textexporter'
include { OPENMS_MZTABEXPORTER as OPENMS_MZTABEXPORTER_QUANT }              from '../../modules/local/openms_mztabexporter'

workflow PROCESS_FEATURE {
    take:
        ch_runs_to_be_quantified

    main:
        ch_versions = Channel.empty()

        // Quantify identifications using targeted feature extraction
        OPENMS_FEATUREFINDERIDENTIFICATION(ch_runs_to_be_quantified)
        ch_versions = ch_versions.mix(OPENMS_FEATUREFINDERIDENTIFICATION.out.versions.first().ifEmpty(null))
        // Link extracted features
        OPENMS_FEATURELINKERUNLABELEDKD(
            OPENMS_FEATUREFINDERIDENTIFICATION.out.featurexml
                .flatMap {
                    meta, raw ->
                        [[[id:meta.sample + "_" + meta.condition, sample:meta.sample, condition:meta.condition, ext:meta.ext], raw]]
                }
                .groupTuple(by:[0]))
        ch_versions = ch_versions.mix(OPENMS_FEATURELINKERUNLABELEDKD.out.versions.first().ifEmpty(null))
        // Resolve conflicting ids matching to the same feature
        OPENMS_IDCONFLICTRESOLVER(OPENMS_FEATURELINKERUNLABELEDKD.out.consensusxml)
        ch_versions = ch_versions.mix(OPENMS_IDCONFLICTRESOLVER.out.versions.first().ifEmpty(null))
        // Export all information as text to csv
        OPENMS_TEXTEXPORTER_QUANT(OPENMS_IDCONFLICTRESOLVER.out.consensusxml)
        ch_versions = ch_versions.mix(OPENMS_TEXTEXPORTER_QUANT.out.versions.first().ifEmpty(null))
        // Export all information as mzTab
        OPENMS_MZTABEXPORTER_QUANT(OPENMS_IDCONFLICTRESOLVER.out.consensusxml)
        ch_versions = ch_versions.mix(OPENMS_MZTABEXPORTER_QUANT.out.versions.first().ifEmpty(null))
    emit:
        // Define the information that is returned by this workflow
        versions = ch_versions
        mztab = OPENMS_MZTABEXPORTER_QUANT.out.mztab
}

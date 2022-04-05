/*
 * Perform the quantification of the samples when the parameter --skip_quantification is not provided
 */

include { OPENMS_FEATUREFINDERIDENTIFICATION }                              from '../../modules/local/openms_featurefinderidentification'
include { OPENMS_FEATURELINKERUNLABELEDKD }                                 from '../../modules/local/openms_featurelinkerunlabeledkd'
include { OPENMS_IDCONFLICTRESOLVER }                                       from '../../modules/local/openms_idconflictresolver'
include { OPENMS_TEXTEXPORTER as OPENMS_TEXTEXPORTER_QUANTIFIED }           from '../../modules/local/openms_textexporter'
include { OPENMS_MZTABEXPORTER }                                            from '../../modules/local/openms_mztabexporter'

workflow POST_QUANTIFICATION {
    take:
        psms_outcome
        aligned_mzml
        filter_q_value

    main:
        ch_versions = Channel.empty()
        // Combining the necessary information into one channel
        psms_outcome
            .join( aligned_mzml, by: [0] )
            .map { it -> [it[0].sample, it[0], it[1], it[2]] }
            .combine( filter_q_value , by: [0] )
            .map { it -> [it[1], it[2], it[3], it[5]] }
            .set{ joined_mzmls_ids_quant }
        // Quantify identifications using targeted feature extraction
        OPENMS_FEATUREFINDERIDENTIFICATION(joined_mzmls_ids_quant)
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
        OPENMS_TEXTEXPORTER_QUANTIFIED(OPENMS_IDCONFLICTRESOLVER.out.consensusxml)
        ch_versions = ch_versions.mix(OPENMS_TEXTEXPORTER_QUANTIFIED.out.versions.first().ifEmpty(null))
        // Export all information as mzTab
        OPENMS_MZTABEXPORTER(OPENMS_IDCONFLICTRESOLVER.out.consensusxml)
        ch_versions = ch_versions.mix(OPENMS_MZTABEXPORTER.out.versions.first().ifEmpty(null))
    emit:
        // Define the information that is returned by this workflow
        versions = ch_versions
        mztab = OPENMS_MZTABEXPORTER.out.mztab
}

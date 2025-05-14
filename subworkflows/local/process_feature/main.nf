/*
 * Perform the quantification by extracting the feature intensities and group runs corresponding to the same sample and condition.
 */

include { OPENMS_FEATUREFINDERIDENTIFICATION }                              from '../../../modules/local/openms/featurefinderidentification'
include { OPENMS_FEATURELINKERUNLABELEDKD }                                 from '../../../modules/local/openmsthirdparty/featurelinkerunlabeledkd'
include { OPENMS_IDCONFLICTRESOLVER }                                       from '../../../modules/local/openms/idconflictresolver'

workflow PROCESS_FEATURE {
    take:
        ch_runs_to_be_quantified

    main:
        ch_versions = Channel.empty()

        // Quantify identifications using targeted feature extraction
        OPENMS_FEATUREFINDERIDENTIFICATION(ch_runs_to_be_quantified).featurexml
                .map { meta, featurexml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), featurexml] }
                .groupTuple()
                .set { ch_features_grouped }
        ch_versions = ch_versions.mix(OPENMS_FEATUREFINDERIDENTIFICATION.out.versions)

        // Link extracted features
        OPENMS_FEATURELINKERUNLABELEDKD(ch_features_grouped)
        ch_versions = ch_versions.mix(OPENMS_FEATURELINKERUNLABELEDKD.out.versions)

        // Resolve conflicting ids matching to the same feature
        OPENMS_IDCONFLICTRESOLVER(OPENMS_FEATURELINKERUNLABELEDKD.out.consensusxml)
        ch_versions = ch_versions.mix(OPENMS_IDCONFLICTRESOLVER.out.versions)

    emit:
        // Define the information that is returned by this workflow
        versions = ch_versions
        consensusxml = OPENMS_IDCONFLICTRESOLVER.out.consensusxml
}

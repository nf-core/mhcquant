/*
 * Perform the quantification by extracting the feature intensities and group runs corresponding to the same sample and condition.
 */

include { OPENMS_FEATUREFINDERIDENTIFICATION }                              from '../../../modules/local/openms/featurefinderidentification'
include { OPENMS_FEATURELINKERUNLABELEDKD }                                 from '../../../modules/local/openmsthirdparty/featurelinkerunlabeledkd'
include { OPENMS_IDCONFLICTRESOLVER }                                       from '../../../modules/local/openms/idconflictresolver/main'
include { OPENMS_FILECONVERTER }                                            from '../../../modules/local/openms/fileconverter/main'

workflow PROCESS_FEATURE {
    take:
        ch_runs_to_be_quantified

    main:
        ch_versions = Channel.empty()

        // Quantify identifications using targeted feature extraction
        OPENMS_FEATUREFINDERIDENTIFICATION(ch_runs_to_be_quantified).featurexml
                .map { meta, featurexml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), featurexml] }
                .groupTuple()
                .set { ch_featuresxmls }
        ch_versions = ch_versions.mix(OPENMS_FEATUREFINDERIDENTIFICATION.out.versions)

        ch_featuresxmls
                .branch {
                    single: it[1].size() == 1
                    multiple: it[1].size() > 1
                }
                .set { ch_features }

        // Link extracted features
        OPENMS_FEATURELINKERUNLABELEDKD(ch_features.multiple)
        ch_versions = ch_versions.mix(OPENMS_FEATURELINKERUNLABELEDKD.out.versions)

        // Single replicate: promote featureXML to consensusXML
        OPENMS_FILECONVERTER(ch_features.single.map { meta, features -> [ meta, features[0], "consensusXML" ] })
        ch_versions = ch_versions.mix(OPENMS_FILECONVERTER.out.versions)

        // Resolve conflicting ids matching to the same feature
        ch_consensus_input = OPENMS_FEATURELINKERUNLABELEDKD.out.consensusxml
            .mix(OPENMS_FILECONVERTER.out.consensusxml)
        
        OPENMS_IDCONFLICTRESOLVER(ch_consensus_input)
        ch_versions = ch_versions.mix(OPENMS_IDCONFLICTRESOLVER.out.versions)

    emit:
        // Define the information that is returned by this workflow
        versions = ch_versions
        consensusxml = OPENMS_IDCONFLICTRESOLVER.out.consensusxml
}

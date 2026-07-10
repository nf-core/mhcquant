/*
 * Perform the quantification by extracting the feature intensities and group runs corresponding to the same sample and condition.
 */

include { OPENMS_FEATUREFINDERIDENTIFICATION } from '../../../modules/local/openms/featurefinderidentification'
include { OPENMS_FEATURELINKERUNLABELEDKD    } from '../../../modules/local/openmsthirdparty/featurelinkerunlabeledkd'
include { OPENMS_IDCONFLICTRESOLVER          } from '../../../modules/local/openms/idconflictresolver/main'
include { OPENMS_FILECONVERTER               } from '../../../modules/nf-core/openms/fileconverter/main'

workflow PROCESS_FEATURE {
    take:
    ch_runs_to_be_quantified

    main:
    // Quantify identifications using targeted feature extraction
    OPENMS_FEATUREFINDERIDENTIFICATION(ch_runs_to_be_quantified).featurexml
        .map { meta, featurexml -> [groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), featurexml] }
        .groupTuple()
        // groupTuple() emits in task-arrival order, which sets the consensus map (rt_N/mz_N/intensity_N)
        // column order. Sort by the leading run ID so per-replicate columns are reproducible across runs.
        .map { key, featurexmls -> [key, featurexmls.sort { it.name.tokenize('_')[0] as int }] }
        .set { ch_featuresxmls }

    ch_featuresxmls
        .branch { _meta, features ->
            single: features.size() == 1
            multiple: features.size() > 1
        }
        .set { ch_features }

    // Link extracted features
    OPENMS_FEATURELINKERUNLABELEDKD(ch_features.multiple)

    // Single replicate: promote featureXML to consensusXML
    OPENMS_FILECONVERTER(ch_features.single.map { meta, features -> [meta, features[0], "consensusXML"] })

    // Resolve conflicting ids matching to the same feature
    ch_consensus_input = OPENMS_FEATURELINKERUNLABELEDKD.out.consensusxml.mix(OPENMS_FILECONVERTER.out.converted)

    OPENMS_IDCONFLICTRESOLVER(ch_consensus_input)

    emit:
    consensusxml = OPENMS_IDCONFLICTRESOLVER.out.consensusxml
}

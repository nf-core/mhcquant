/*
 * Perform an additional step where the process are collected
 * that are called when the parameter "refine_fdr" is provided
 */

include { OPENMS_MZTABEXPORTER as OPENMS_MZTABEXPORTERPERC } from '../../modules/local/openms_mztabexporter'
include { OPENMS_MZTABEXPORTER as OPENMS_MZTABEXPORTERPSM }  from '../../modules/local/openms_mztabexporter'
include { MHCFLURRY_PREDICTPSMS }                            from '../../modules/local/mhcflurry_predictpsms'
include { OPENMS_PERCOLATORADAPTER }                         from '../../modules/local/openms_percolatoradapter'
include { OPENMS_IDFILTER as OPENMS_IDFILTER_PSMS }          from '../../modules/local/openms_idfilter'
include { OPENMS_IDFILTER as OPENMS_IDFILTER_REFINED }       from '../../modules/local/openms_idfilter'

workflow REFINE_FDR {
    // Define the input parameters
    take:
        filtered_perc_output
        psm_features
        classI_alleles

    main:
        ch_versions = Channel.empty()
        // Export filtered percolator results as mztab
        OPENMS_MZTABEXPORTERPERC( filtered_perc_output )
        ch_versions = ch_versions.mix(OPENMS_MZTABEXPORTERPERC.out.versions)
        // Export psm results as mztab
        OPENMS_MZTABEXPORTERPSM( psm_features )
        ch_versions = ch_versions.mix(OPENMS_MZTABEXPORTERPSM.out.versions)
        // Predict psm results using mhcflurry to shrink search space
        MHCFLURRY_PREDICTPSMS(
            OPENMS_MZTABEXPORTERPERC.out.mztab
                .join( OPENMS_MZTABEXPORTERPSM.out.mztab, by:[0] )
                .map{ it -> [it[0].sample, it[0], it[1], it[2]] }
                .combine( classI_alleles, by:0)
                .map(it -> [it[1], it[2], it[3], it[4]])
        )
        ch_versions = ch_versions.mix(MHCFLURRY_PREDICTPSMS.out.versions)

        // Filter psm results by shrinked search space
        OPENMS_IDFILTER_PSMS(psm_features.combine( MHCFLURRY_PREDICTPSMS.out.idxml, by: [0] ))
        ch_versions = ch_versions.mix(OPENMS_IDFILTER_PSMS.out.versions)
        // Recompute percolator fdr on shrinked search space
        OPENMS_PERCOLATORADAPTER( OPENMS_IDFILTER_PSMS.out.idxml )
        ch_versions = ch_versions.mix(OPENMS_PERCOLATORADAPTER.out.versions)
        // Filter results by refined fdr
        OPENMS_IDFILTER_REFINED(OPENMS_PERCOLATORADAPTER.out.idxml.flatMap { it -> [tuple(it[0], it[1], null)]})
        ch_versions = ch_versions.mix(OPENMS_IDFILTER_REFINED.out.versions)
    emit:
        // Define the information that is returned by this workflow
        filter_refined_q_value = OPENMS_IDFILTER_REFINED.out.idxml
        versions = ch_versions
}

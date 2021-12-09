/*
 * Perform an additional step where the process are collected
 * that are called when the paramater "refine_fdr_on_predicted_subset" is provided
 */

params.exporter_prec_options  = [:]
params.exporter_psm_options  = [:]
params.percolator_adapter_refine_options = [:]
params.whitelist_filter_options = [:]
params.filter_options = [:]

def openms_mztab_exporter_prec_options = params.exporter_prec_options.clone()
def openms_mztab_exporter_psm_options = params.exporter_psm_options.clone()
def openms_percolator_adapter_options = params.percolator_adapter_refine_options.clone()
def openms_id_filter_psms_options = params.whitelist_filter_options.clone()
def openms_id_filter_qvalue_options = params.filter_options.clone()

include { OPENMS_MZTABEXPORTER as OPENMS_MZTABEXPORTERPERC } from '../../modules/local/openms_mztabexporter'                                       addParams( options: openms_mztab_exporter_prec_options )
include { OPENMS_MZTABEXPORTER as OPENMS_MZTABEXPORTERPSM }  from '../../modules/local/openms_mztabexporter'                                       addParams( options: openms_mztab_exporter_psm_options )
include { MHCFLURRY_PREDICTPSMS }                            from '../../modules/local/mhcflurry_predictpsms'                                      addParams( options: [:] )
include { OPENMS_PERCOLATORADAPTER }                         from '../../modules/local/openms_percolatoradapter'                                   addParams( options: openms_percolator_adapter_options )
include { OPENMS_IDFILTER as OPENMS_IDFILTER_PSMS }          from '../../modules/local/openms_idfilter'                                            addParams( options: openms_id_filter_psms_options )
include { OPENMS_IDFILTER as OPENMS_IDFILTER_REFINED }       from '../../modules/local/openms_idfilter'                                            addParams( options: openms_id_filter_qvalue_options )

workflow REFINE_FDR_ON_PREDICTED_SUBSET {
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

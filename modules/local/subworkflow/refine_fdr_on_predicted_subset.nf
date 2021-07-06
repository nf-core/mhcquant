/*
 * Perform an additional step where the process are collected 
 * that are called when the paramater "refine_fdr_on_predicted_subset" is provided
 */

// VALIDATED (EQUAL TO THE OLD CODE)

params.percolator_adapter_options = [:]
params.filter_options = [:]
params.whitelist_filter_options = [:]

def percolator_adapter_options = params.percolator_adapter_options.clone()
percolator_adapter_options.suffix = "perc_subset"

def filter_psms_options = params.whitelist_filter_options.clone()
def filter_refined_qvalue_options = params.filter_options.clone()

filter_psms_options.suffix = "pred_filtered"
filter_refined_qvalue_options.suffix = "perc_subset_filtered"

include { OPENMS_MZTABEXPORTER as OPENMS_MZTABEXPORTERPERC } from '../openms_mztabexporter'                                       addParams( options: [ suffix: "all_ids_merged_psm_perc_filtered" ] )
include { OPENMS_MZTABEXPORTER as OPENMS_MZTABEXPORTERPSM }  from '../openms_mztabexporter'                                       addParams( options: [ suffix: "all_ids_merged" ] )
include { PREDICT_PSMS }                                     from '../predict_psms'                                               addParams( options: [:] )
include { OPENMS_PERCOLATORADAPTER }                         from '../openms_percolatoradapter'                                   addParams( options: percolator_adapter_options )
include { OPENMS_IDFILTER as OPENMS_IDFILTER_PSMS }          from '../openms_idfilter'                                            addParams( options: filter_psms_options )
include { OPENMS_IDFILTER as OPENMS_IDFILTER_REFINED }       from '../openms_idfilter'                                            addParams( options: filter_refined_qvalue_options )



workflow REFINE_FDR_ON_PREDICTED_SUBSET {
    // Define the input parameters
    take:
        filtered_perc_output
        psm_features
        classI_alleles

    main:
        ch_software_versions = Channel.empty()
        // Export filtered percolator results as mztab
        OPENMS_MZTABEXPORTERPERC( filtered_perc_output )
        ch_software_versions = ch_software_versions.mix(OPENMS_MZTABEXPORTERPERC.out.version.first().ifEmpty(null))
        // Export psm results as mztab
        OPENMS_MZTABEXPORTERPSM( psm_features )
        // Predict psm results using mhcflurry to shrink search space
        PREDICT_PSMS(OPENMS_MZTABEXPORTERPERC.out.mztab.join( OPENMS_MZTABEXPORTERPSM.out.mztab, by:[0,1] ).combine( classI_alleles, by:1 ) )
        // Filter psm results by shrinked search space
        OPENMS_IDFILTER_PSMS(psm_features.combine( PREDICT_PSMS.out.idxml, by: [0, 1] ))
        // Replace the id's and add the condition (place fillers)
        ch_predicted_psms = OPENMS_IDFILTER_PSMS.out.idxml.flatMap { it -> [ tuple( "id", it[1], "condition", it[3] ) ] }
        // Recompute percolator fdr on shrinked search space
        OPENMS_PERCOLATORADAPTER( ch_predicted_psms ) 
        // Filter results by refined fdr
        OPENMS_IDFILTER_REFINED(OPENMS_PERCOLATORADAPTER.out.idxml.flatMap { it -> [tuple(it[0], it[1], it[2], it[3], null)]})
        
    emit:
        // Define the information that is returned by this workflow
        filter_refined_q_value = OPENMS_IDFILTER_REFINED.out.idxml
        version = ch_software_versions
}
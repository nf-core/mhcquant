/*
 * Perform an additional step where the process are collected 
 * that are called when the paramater "refine_fdr_on_predicted_subset" is provided
 */

// VALIDATED (EQUAL TO THE OLD CODE)

params.percolator_adapter_options = [:]
def percolator_adapter_options = params.percolator_adapter_options.clone()
percolator_adapter_options.suffix = "perc_subset"

include { OPENMS_MZTABEXPORTER as OPENMS_MZTABEXPORTERPERC } from '../process/openms_mztabexporter'                                       addParams( options: [ suffix: "all_ids_merged_psm_perc_filtered" ] )
include { OPENMS_MZTABEXPORTER as OPENMS_MZTABEXPORTERPSM }  from '../process/openms_mztabexporter'                                       addParams( options: [ suffix: "all_ids_merged" ] )

include { PREDICT_PSMS }                                     from '../process/predict_psms'                                               addParams( options: [:] )
include { FILTER_PSMS_BY_PREDICTIONS }                       from '../process/filter_psms_by_predictions'                                 addParams( options: [:] )
include { OPENMS_PERCOLATORADAPTER }                         from '../process/openms_percolatoradapter'                                  addParams( options: percolator_adapter_options )
include { FILTER_REFINED_Q_VALUE }                           from '../process/filter_refined_q_value'                                     addParams( options: [:] )

workflow REFINE_FDR_ON_PREDICTED_SUBSET {
    // Define the input parameters
    take:
        filtered_perc_output
        psm_features
        classI_alleles

    main:
        // Export filtered percolator results as mztab
        OPENMS_MZTABEXPORTERPERC( filtered_perc_output ) // Include an if to subtitude for the when: params.refine_fdr_on_predicted_subset
        // Export psm results as mztab
        OPENMS_MZTABEXPORTERPSM( psm_features )
        // Predict psm results using mhcflurry to shrink search space
        // .combine( classI_alleles, by:1 ).view()
        PREDICT_PSMS(OPENMS_MZTABEXPORTERPERC.out.mztab.join( OPENMS_MZTABEXPORTERPSM.out.mztab, by:[0,1] ).combine( classI_alleles, by:1 ) )
        // Filter psm results by shrinked search space
        FILTER_PSMS_BY_PREDICTIONS( psm_features, PREDICT_PSMS.out.idxml )
        // Replace the id's and add the condition (place fillers)
        ch_predicted_psms = FILTER_PSMS_BY_PREDICTIONS.out.idxml.flatMap { it -> [ tuple( "id", it[1], "condition", it[2] ) ] }
        // Recompute percolator fdr on shrinked search space
        OPENMS_PERCOLATORADAPTER( ch_predicted_psms ) 
        // Filter results by refined fdr
        FILTER_REFINED_Q_VALUE( OPENMS_PERCOLATORADAPTER.out.idxml.flatMap { it -> [ tuple( it[0], it[1], it[3] ) ] })

    emit:
        // Define the information that is returned by this workflow
        filter_refined_q_value = FILTER_REFINED_Q_VALUE.out.idxml
}
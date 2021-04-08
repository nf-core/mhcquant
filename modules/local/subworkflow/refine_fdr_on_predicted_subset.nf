/*
 * Perform an additional step where the process are collected 
 * that are called when the paramater "refine_fdr_on_predicted_subset" is provided
 */

params.options = [:]

// include { EXPORT_MZTAB_PERC }                               from '../process/export_mztab_perc'                                          addParams( options: [:] )
// include { EXPORT_MZTAB_PSM }                                from '../process/export_mztab_psm'                                           addParams( options: [:] )
include { OPENMS_MZTABEXPORTER as OPENMS_MZTABEXPORTERPERC } from '../process/openms_mztabexporter'                                       addParams( options: [ suffix: "all_ids_merged_psm_perc_filtered" ] )
include { OPENMS_MZTABEXPORTER as OPENMS_MZTABEXPORTERPSM }  from '../process/openms_mztabexporter'                                       addParams( options: [ suffix: "all_ids_merged" ] )


include { PREDICT_PSMS }                                     from '../process/predict_psms'                                               addParams( options: [:] )
include { FILTER_PSMS_BY_PREDICTIONS }                       from '../process/filter_psms_by_predictions'                                 addParams( options: [:] )
include { RUN_PERCOLATOR_ON_PREDICTED_SUBSET }               from '../process/run_percolator_on_predicted_subset'                         addParams( options: [:] )
include { FILTER_REFINED_Q_VALUE }                           from '../process/filter_refined_q_value'                                     addParams( options: [:] )

workflow REFINE_FDR_ON_PREDICTED_SUBSET {
    // Define the input parameters
    take:
        filtered_perc_output
        psm_features
        classI_alleles
        fdr_level // find another way

    main:
        // Export filtered percolator results as mztab
        OPENMS_MZTABEXPORTERPERC( filtered_perc_output ) // Include an if to subtitude for the when: params.refine_fdr_on_predicted_subset
        // Export psm results as mztab
        // psm_features.view()
        // OPENMS_MZTABEXPORTERPERC.out.mztab.view()
        OPENMS_MZTABEXPORTERPSM( psm_features )
        // OPENMS_MZTABEXPORTERPSM.out[0].view()
        OPENMS_MZTABEXPORTERPERC.out.mztab.join( OPENMS_MZTABEXPORTERPSM.out.mztab, by:[0,1] ).combine( classI_alleles, by:1 ).view()
        // Predict psm results using mhcflurry to shrink search space
        // TODO: Create a new container with an older python version 2.7.15
        PREDICT_PSMS(OPENMS_MZTABEXPORTERPERC.out.mztab.join( OPENMS_MZTABEXPORTERPSM.out.mztab, by:[0,1] ).combine( classI_alleles, by:1 ) )
        // Filter psm results by shrinked search space
        FILTER_PSMS_BY_PREDICTIONS( psm_features, PREDICT_PSMS.out.idxml )
        // Recompute percolator fdr on shrinked search space
        RUN_PERCOLATOR_ON_PREDICTED_SUBSET( FILTER_PSMS_BY_PREDICTIONS.out.idxml, fdr_level ) 
        // Filter results by refined fdr
        FILTER_REFINED_Q_VALUE( RUN_PERCOLATOR_ON_PREDICTED_SUBSET.out.idxml )

    emit:
        // Define the information that is returned by this workflow
        filter_refined_q_value = FILTER_REFINED_Q_VALUE.out.idxml
}
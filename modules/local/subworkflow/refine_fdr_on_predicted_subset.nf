/*
 * Perform an additional step where the process are collected 
 * that are called when the paramater "refine_fdr_on_predicted_subset" is provided
 */

params.options = [:]

include { EXPORT_MZTAB_PERC }                               from '../process/export_mztab_perc'                                          addParams( options: [:] )
include { EXPORT_MZTAB_PSM }                                from '../process/export_mztab_psm'                                           addParams( options: [:] )
include { PREDICT_PSMS }                                    from '../process/predict_psms'                                               addParams( options: [:] )
include { FILTER_PSMS_BY_PREDICTIONS }                      from '../process/filter_psms_by_predictions'                                 addParams( options: [:] )
include { RUN_PERCOLATOR_ON_PREDICTED_SUBSET }              from '../process/run_percolator_on_predicted_subset'                         addParams( options: [:] )
include { FILTER_REFINED_Q_VALUE }                          from '../process/filter_refined_q_value'                                     addParams( options: [:] )

workflow REFINE_FDR_ON_PREDICTED_SUBSET {
    // Define the input parameters
    take:
        filtered_perc_output
        psm_features
        classI_alleles
        fdr_level
        filtered_ids 
        aligned_mzml_files

    main:
        // Export filtered percolator results as mztab
        EXPORT_MZTAB_PERC(filtered_perc_output) 
        // Export psm results as mztab
        EXPORT_MZTAB_PSM(psm_features)
        // Predict psm results using mhcflurry to shrink search space
        PREDICT_PSMS(EXPORT_MZTAB_PERC.out.join(EXPORT_MZTAB_PSM.out, by:[0,1]).combine(classI_alleles, by:1) )
        // Filter psm results by shrinked search space
        FILTER_PSMS_BY_PREDICTIONS(psm_features, PREDICT_PSMS.out)
        // Recompute percolator fdr on shrinked search space
        RUN_PERCOLATOR_ON_PREDICTED_SUBSET(FILTER_PSMS_BY_PREDICTIONS.out, fdr_level) 
        // Filter results by refined fdr
        FILTER_REFINED_Q_VALUE(RUN_PERCOLATOR_ON_PREDICTED_SUBSET.out)

        filtered_ids
           .flatMap { it -> [tuple(it[0], it[1], it[2], it[3])]}
           .join(aligned_mzml_files, by: [0,1,2])
           .combine(filtered_perc_output.mix(FILTER_REFINED_Q_VALUE.out), by:1)
           .set{joined_mzmls_ids_quant}

    emit:
        // Define the information that is returned by this workflow
        joined_mzmls_ids_quant
}
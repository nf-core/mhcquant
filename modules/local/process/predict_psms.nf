// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

// TODO: add python 2.7.15 and mhcflurry 1.4.3
process PREDICT_PSMS {
    publishDir "${params.outdir}/Intermediate_Results/"

    input:
    tuple val(Sample), val(id), file(perc_mztab_file), file(psm_mztab_file), val(d), val(allotypes_refine)

    output:
    tuple val("$id"), val("$Sample"), file("${Sample}_peptide_filter.idXML")

    when:
    params.refine_fdr_on_predicted_subset

    script:
    """
        mhcflurry-downloads --quiet fetch models_class1
        mhcflurry_predict_mztab_for_filtering.py ${params.subset_affinity_threshold} '${allotypes_refine}' ${perc_mztab_file} ${psm_mztab_file} ${Sample}_peptide_filter.idXML
    """
}

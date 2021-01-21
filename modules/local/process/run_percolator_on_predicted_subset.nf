// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process RUN_PERCOLATOR_ON_PREDICTED_SUBSET {
    publishDir "${params.outdir}/Intermediate_Results/"

    input:
    tuple val(id), val(Sample), file(id_file_psm_subset)
    val fdr_level

    output:
    tuple val("$id"), val("$Sample"), file("${Sample}_perc_subset.idXML")

    when:
    params.refine_fdr_on_predicted_subset

    script:
    """
        OMP_NUM_THREADS=${task.cpus} \\
        PercolatorAdapter -in ${id_file_psm_subset} \\
            -out ${Sample}_perc_subset.idXML \\
            -seed 4711 \\
            -trainFDR 0.05 \\
            -testFDR 0.05 \\
            -enzyme no_enzyme \\
            -subset-max-train ${params.subset_max_train} \\
            -doc ${params.description_correct_features} \\
            $fdr_level
    """
}
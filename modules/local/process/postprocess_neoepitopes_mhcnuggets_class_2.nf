// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process POSTPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2 {
    publishDir "${params.outdir}/class_2_bindings"

    input:
    tuple val(Sample), val(id), file(neoepitopes), val(d), file(predicted_cl_2) 

    output:
    tuple val("$Sample"), file("*.csv")

    when:
    params.include_proteins_from_vcf & params.predict_class_2

    script:
    """
        postprocess_neoepitopes_mhcnuggets.py --input ${predicted_cl_2} --neoepitopes ${neoepitopes}
    """
}
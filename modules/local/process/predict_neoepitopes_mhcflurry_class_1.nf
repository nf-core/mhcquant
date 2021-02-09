// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

// TODO: add python 2.7.15 and mhcflurry 1.4.3
process PREDICT_NEOEPITOPES_MHCFLURRY_CLASS_1 {
    publishDir "${params.outdir}/class_1_bindings"

    input:
    tuple val(Sample), val(id), val(allotypes), val(d), file(neoepitopes) 

    output:
    tuple val("$id"), val("$Sample"), file("*_${Sample}_predicted_neoepitopes_class_1.csv")
    
    when:
    params.include_proteins_from_vcf & params.predict_class_1

    script:
    """
        mhcflurry-downloads --quiet fetch models_class1
        mhcflurry_neoepitope_binding_prediction.py '${allotypes}' ${neoepitopes} _${Sample}_predicted_neoepitopes_class_1.csv
    """
}
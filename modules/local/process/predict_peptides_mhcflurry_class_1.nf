// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

// TODO: add python 2.7.15 and mhcflurry 1.4.3
process PREDICT_PEPTIDES_MHCFLURRY_CLASS_1 {
    publishDir "${params.outdir}/class_1_bindings"

    input:
    tuple val(Sample), val(id), file(mztab_file), val(d), val(class_1_alleles)

    output:
    tuple val("$Sample"), file("*predicted_peptides_class_1.csv")

    when:
    params.predict_class_1

    script:
    """
        mhcflurry-downloads --quiet fetch models_class1
        mhcflurry_predict_mztab.py '${class_1_alleles}' ${mztab_file} ${Sample}_predicted_peptides_class_1.csv
    """
}
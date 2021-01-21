// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process PREDICT_NEOEPITOPES_MHCNUGGETS_CLASS_2 {

    input:
    tuple val(Sample), val(id), file(preprocessed_neoepitopes), val(d), val(cl_2_alleles)

    output:
    tuple val("$id"), val("$Sample"), file("*_predicted_neoepitopes_class_2")

    when:
    params.include_proteins_from_vcf & params.predict_class_2

    script:
    """
        mhcnuggets_predict_peptides.py --peptides ${preprocessed_neoepitopes} --alleles '${cl_2_alleles}' --output _predicted_neoepitopes_class_2
    """
}
// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2 {

    input:
    tuple val(Sample), val(id), file(preprocessed_peptides), val(d), val(class_2_alleles)

    output:
    tuple val("$id"), val("$Sample"), file("*_predicted_peptides_class_2")

    when:
    params.predict_class_2

    script:
    """
    mhcnuggets_predict_peptides.py --peptides ${preprocessed_peptides} --alleles '${class_2_alleles}' --output _predicted_peptides_class_2
    """
}
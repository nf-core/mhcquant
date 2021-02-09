// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

// TODO: add python 2.7.15 and mhcnuggets 2.3.2
process PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2 {

    input:
    tuple val(id), val(Sample), file(neoepitopes)

    output:
    tuple val("$id"), val("$Sample"), file("${Sample}_mhcnuggets_preprocessed")

    when:
    params.include_proteins_from_vcf & params.predict_class_2

    script:
    """
    preprocess_neoepitopes_mhcnuggets.py --neoepitopes ${neoepitopes} --output ${Sample}_mhcnuggets_preprocessed
    """
}

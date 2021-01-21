// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2 {

    input:
    tuple val(id), val(Sample), file(mztab_file) 

    output:
    tuple val("$id"), val("$Sample"), file("${Sample}_preprocessed_mhcnuggets_peptides")
    tuple val("$id"), val("$Sample"), file('peptide_to_geneID') 

    // emit:
    // preprocessed_mhcnuggets_peptides
    // peptide_to_geneID
// 
    when:
    params.predict_class_2

    script:
    """
        preprocess_peptides_mhcnuggets.py --mztab ${mztab_file} --output ${Sample}_preprocessed_mhcnuggets_peptides
    """
}

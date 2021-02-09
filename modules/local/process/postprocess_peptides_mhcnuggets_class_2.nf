// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

// TODO: add python 2.7.15 and mhcnuggets 2.3.2
process POSTPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2 {
    publishDir "${params.outdir}/class_2_bindings"

    input:
    tuple val(Sample), val(id), file(predicted_peptides), val(d), file(peptide_to_geneID)

    output:
    tuple val("$Sample"), file('*.csv')

    when:
    params.predict_class_2

    script:
    """
        postprocess_peptides_mhcnuggets.py --input ${predicted_peptides} --peptides_seq_ID ${peptide_to_geneID} --output ${Sample}_postprocessed.csv
    """
}

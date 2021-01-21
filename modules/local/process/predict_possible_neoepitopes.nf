// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process PREDICT_POSSIBLE_NEOEPITOPES {
    publishDir "${params.outdir}/"

    label 'process_web'

    input:
    tuple val(id), val(Sample), val(alleles), file(vcf_file)

    output:
    tuple val("id"), val("$Sample"), file("${Sample}_vcf_neoepitopes.csv")
    tuple val("id"), val("$Sample"), file("${Sample}_vcf_neoepitopes.txt")

    when:
    params.include_proteins_from_vcf & params.predict_class_1

    script:
    """
        vcf_neoepitope_predictor.py -t ${params.variant_annotation_style} -r ${params.variant_reference} -a '${alleles}' -minl ${params.peptide_min_length} -maxl ${params.peptide_max_length} -v ${vcf_file} -o ${Sample}_vcf_neoepitopes.csv
    """
}
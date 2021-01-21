// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process GENERATE_PROTEINS_FROM_VCF {
    publishDir "${params.outdir}/"

    label 'process_web'

    input:
    tuple val(Sample), val(id), file(fasta_file_vcf), val(d), file(vcf_file)
    
    output:
    tuple val("$id"), val("$Sample"), file("${Sample}_${fasta_file_vcf.baseName}_added_vcf.fasta")

    when:
    params.include_proteins_from_vcf

    script:
    """
    variants2fasta.py -v ${vcf_file} -t ${params.variant_annotation_style} -r ${params.variant_reference} -f ${fasta_file_vcf} -o ${Sample}_${fasta_file_vcf.baseName}_added_vcf.fasta ${variant_indel_filter} ${variant_snp_filter} ${variant_frameshift_filter}
    """
}
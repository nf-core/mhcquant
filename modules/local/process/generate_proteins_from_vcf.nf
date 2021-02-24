// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

def VERSION.FRED2 = '2.0.6'
def VERSION.MHCNUGGETS = '2.3.2'

//TODO: combine in a subflow --> when needs to be removed
process GENERATE_PROTEINS_FROM_VCF {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'.', publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    }

    input:
        tuple val(Sample), val(id), file(fasta_file), val(d), file(vcf_file)
        val variant_indel_filter
        val variant_snp_filter
        val variant_frameshift_filter
    
    output:
        tuple val("$id"), val("$Sample"), file("${Sample}_${fasta_file.baseName}_added_vcf.fasta"), emit: vcf_fasta
        path  "*.version.txt", emit: version

    when:
        params.include_proteins_from_vcf

    script:
    """
        variants2fasta.py -v ${vcf_file} -t ${params.variant_annotation_style} -r ${params.variant_reference} -f ${fasta_file} -o ${Sample}_${fasta_file.baseName}_added_vcf.fasta ${variant_indel_filter} ${variant_snp_filter} ${variant_frameshift_filter}
        
        echo $VERSION.FRED2 > fred2.version.txt
        echo $VERSION.MHCNUGGETS > mhcnuggets.version.txt
        mhcflurry-predict --version &> mhcflurry.version.txt
    """
}
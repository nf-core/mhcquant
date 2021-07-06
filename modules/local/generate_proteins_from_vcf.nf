// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]
options    = initOptions(params.options)

def VERSIONFRED2 = '2.0.6'
def VERSIONMHCNUGGETS = '2.3.2'

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
        tuple val(Sample), val(id), path(fasta_file), val(d), path(vcf_file)
    
    output:
        tuple val("$id"), val("$Sample"), path("*_vcf.fasta"), emit: vcf_fasta
        path  "*.version.txt", emit: version

    when:
        params.include_proteins_from_vcf

    script:
    """
        variants2fasta.py -v ${vcf_file} -f ${fasta_file} -o ${Sample}_${fasta_file.baseName}_added_vcf.fasta $options.args
        
        echo $VERSIONFRED2 > fred2.version.txt
        echo $VERSIONMHCNUGGETS > mhcnuggets.version.txt
        echo \$(mhcflurry-predict --version 2>&1) | sed 's/^.*mhcflurry //; s/ .*\$//' &> mhcflurry.version.txt
    """
}
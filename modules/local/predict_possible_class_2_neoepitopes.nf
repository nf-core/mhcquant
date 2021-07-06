// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

def VERSIONFRED2 = '2.0.6'
def VERSIONMHCNUGGETS = '2.3.2'

process PREDICT_POSSIBLE_CLASS_2_NEOEPITOPES {
    label 'process_web'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'.', publish_id:'') }

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    }

    input:
        tuple val(id), val(Sample), val(alleles_II), path(vcf_file)

    output:
        tuple val("id"), val("$Sample"), path("${Sample}_vcf_neoepitopes.csv"), emit: csv
        tuple val("id"), val("$Sample"), path("${Sample}_vcf_neoepitopes.txt"), emit: txt   
        path  "*.version.txt", emit: version

    script:
        """
            vcf_neoepitope_predictor.py -t ${params.variant_annotation_style} -r ${params.variant_reference} -a '${alleles_II}' -minl ${params.peptide_min_length} -maxl ${params.peptide_max_length} -v ${vcf_file} -o ${Sample}_vcf_neoepitopes.csv
            echo $VERSIONFRED2 > fred2.version.txt
            echo $VERSIONMHCNUGGETS > mhcnuggets.version.txt
            mhcflurry-predict --version &> mhcflurry.version.txt
        """
}
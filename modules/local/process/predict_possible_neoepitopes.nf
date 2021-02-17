// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

def VERSION = '2.0.7'
// TODO: add python 2.7.15 and fred 2.0.6
//TODO: combine in a subflow --> when needs to be removed
process PREDICT_POSSIBLE_NEOEPITOPES {
    label 'process_web'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'.', publish_id:'') }

    conda (params.enable_conda ? "bioconda::fred2=2.0.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fred2:2.0.7--py_0"
    } else {
        container "quay.io/biocontainers/fred2:2.0.7--py_0"
    }

    input:
        tuple val(id), val(Sample), val(alleles), file(vcf_file)

    output:
        tuple val("id"), val("$Sample"), file("${Sample}_vcf_neoepitopes.csv"), emit: csv   
        tuple val("id"), val("$Sample"), file("${Sample}_vcf_neoepitopes.txt"), emit: txt   
        path  "*.version.txt", emit: version

    when:
        params.include_proteins_from_vcf & params.predict_class_1

    script:
        """
            vcf_neoepitope_predictor.py -t ${params.variant_annotation_style} -r ${params.variant_reference} -a '${alleles}' -minl ${params.peptide_min_length} -maxl ${params.peptide_max_length} -v ${vcf_file} -o ${Sample}_vcf_neoepitopes.ch_software_versions
            
            echo $VERSION > fred2.version.txt
        """
}
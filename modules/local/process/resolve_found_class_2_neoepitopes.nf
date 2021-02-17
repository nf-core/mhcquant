// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

def VERSION = '2.0.7'
//TODO: combine in a subflow --> when needs to be removed
process RESOLVE_FOUND_CLASS_2_NEOEPITOPES {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'.', publish_id:'') }

    echo true

    conda (params.enable_conda ? "bioconda::fred2=2.0.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fred2:2.0.7--py_0"
    } else {
        container "quay.io/biocontainers/fred2:2.0.7--py_0"
    }

    input:
        tuple val(id), val(Sample), file(mztab), file(neoepitopes)

    output:
        tuple val("$id"), val("$Sample"), file("${Sample}_found_neoepitopes_class_2.csv"), emit: csv   
        path  "*.version.txt", emit: version

    when:
        params.include_proteins_from_vcf & params.predict_class_2

    script:
        """
            resolve_neoepitopes.py -n ${neoepitopes} -m ${mztab} -f csv -o ${Sample}_found_neoepitopes_class_2

            echo $VERSION > fred2.version.txt
        """
}
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

def VERSIONFRED2 = '2.0.6'
def VERSIONMHCNUGGETS = '2.3.2'

process RESOLVE_FOUND_CLASS_2_NEOEPITOPES {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'.', publish_id:'') }

    echo true

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    }

    input:
        tuple val(id), val(Sample), path(mztab), path(neoepitopes)

    output:
        tuple val("$id"), val("$Sample"), path("${Sample}_found_neoepitopes_class_2.csv"), emit: csv   
        path  "*.version.txt", emit: version

    script:
        """
            resolve_neoepitopes.py -n ${neoepitopes} -m ${mztab} -f csv -o ${Sample}_found_neoepitopes_class_2
            echo $VERSIONFRED2 > fred2.version.txt
            echo $VERSIONMHCNUGGETS > mhcnuggets.version.txt
            mhcflurry-predict --version &> mhcflurry.version.txt
        """
}
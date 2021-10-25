// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSIONFRED2 = '2.0.6'
def VERSIONMHCNUGGETS = '2.3.2'

process PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2 {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    }

    input:
        tuple val(meta), path(peptides), val(alleles)

    output:
        tuple val(meta), path("*_predicted_peptides_class_2"), emit: csv
        path "versions.yml", emit: versions

    script:
        def prefix = options.suffix ? "${meta.sample}_${options.suffix}" : "${meta.sample}_predicted_peptides_class_2"

        """
            mhcnuggets_predict_peptides.py --peptides ${peptides} --alleles '${alleles}' --output ${prefix}

            cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
                mhcnuggets: \$(echo $VERSIONMHCNUGGETS)
                FRED2: \$(echo $VERSIONFRED2)
            END_VERSIONS
        """
}

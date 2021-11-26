// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MHCNUGGETS_PEPTIDESCLASS2PRE {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0"
    } else {
        container "quay.io/biocontainers/mhcnuggets:2.3.2--py_0"
    }

    input:
        tuple val(meta), path(mztab)

    output:
        tuple val(meta), path("${prefix}")          , emit: preprocessed
        tuple val(meta), path('*peptide_to_geneID*'), emit: geneID
        path "versions.yml"                         , emit: versions

    script:
        def prefix = options.suffix ? "${meta.sample}_${options.suffix}" : "${meta.sample}_preprocessed_mhcnuggets_peptides"

        """
        preprocess_peptides_mhcnuggets.py --mztab $mztab \\
            --output ${prefix}
            
        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
        END_VERSIONS
        """
}

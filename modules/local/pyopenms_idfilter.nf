process PYOPENMS_IDFILTER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pyopenms=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyopenms:3.0.0--py311h9b8898c_0' :
        'biocontainers/pyopenms:3.0.0--py311h9b8898c_0' }"

    input:
        tuple val(meta), path(idxml), path(whitelist)

    output:
        tuple val(meta), path("*_fdr_filtered.idXML") , emit: filtered
        path "versions.yml"                           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        prefix =  task.ext.prefix ?: "${meta.id}_${meta.sample}_${meta.condition}_fdr_filtered"

        """
        IDFilter.py \\
            --input $idxml \\
            --whitelist $whitelist \\
            --output ${prefix}.idXML

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pyopenms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}

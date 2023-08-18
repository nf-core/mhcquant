process PYOPENMS_FDRFILTERRUNS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::pyopenms=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyopenms:3.0.0--py311h9b8898c_0' :
        'biocontainers/pyopenms:3.0.0--py311h9b8898c_0' }"

    input:
        tuple val(merge_id), val(meta), path(idxml), path(pout_idxml)

    output:
        tuple val(merge_id), val(meta), path("*_fdr_filtered.idXML") , emit: fdr_filtered_runs
        path "versions.yml"                           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        prefix =  task.ext.prefix ?: "${meta.id}_${merge_id.id}_fdr_filtered"

        """
        fdr_filter_runs.py \\
            --input $idxml \\
            --pout $pout_idxml \\
            --output ${prefix}.idXML \\

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pyopenms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}

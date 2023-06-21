process OPENMS_PSMFEATUREEXTRACTOR {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::openms=2.8.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.9.1--h135471a_1' :
        'biocontainers/openms:2.9.1--h135471a_1' }"

    input:
        tuple val(meta), path(merged)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${merged.baseName}_psm"
        def args             = task.ext.args ?: ''
        def extra_features = ""
        if(params.use_deeplc){
            def extra_features = "-extra"
            if(params.add_abs_rt_error){
                extra_features = "${extra_features} deeplc_abs_error"
            }
            if(params.add_log_rt_error){
                extra_features = "${extra_features} deeplc_log_error"
            }
            if(params.add_sqr_rt_error || (!params.add_sqr_rt_error && !params.add_abs_rt_error && !params.add_log_rt_error)){
                extra_features = "${extra_features} deeplc_sqr_error"
            }
        }

        """
        PSMFeatureExtractor -in $merged \\
            -out ${prefix}.idXML \\
            -threads $task.cpus \\
            $extra_features \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}

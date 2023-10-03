process OPENMS_PSMFEATUREEXTRACTOR {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::openms=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.0.0--h8964181_1' :
        'biocontainers/openms:3.0.0--h8964181_1' }"

    input:
        tuple val(meta), path(idxml)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${meta.id}_psm"
        def args             = task.ext.args ?: ''
        def extra_features = ""
        if(params.use_deeplc || params.use_ms2pip){
            extra_features = "-extra"
        }
        if(params.use_deeplc){
            if(params.deeplc_add_abs_rt_error){
                extra_features = "${extra_features} deeplc_abs_error"
            }
            if(params.deeplc_add_log_rt_error){
                extra_features = "${extra_features} deeplc_log_error"
            }
            if(params.deeplc_add_sqr_rt_error || (!params.deeplc_add_sqr_rt_error && !params.deeplc_add_abs_rt_error && !params.deeplc_add_log_rt_error)){
                extra_features = "${extra_features} deeplc_sqr_error"
            }
        }
        if(params.use_ms2pip){
            extra_features = "${extra_features} spectrum_correlation"
        }

        """
        PSMFeatureExtractor -in $idxml \\
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

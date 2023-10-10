process MS2PIP {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ms2pip=3.11.0 bioconda::pyopenms=2.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-beb85d5ee68ba9251d26079ca28797d51ea3c49a:857e5e7908422b6ea5016a3c313f67087fbe2f8b-0' :
        'biocontainers/mulled-v2-beb85d5ee68ba9251d26079ca28797d51ea3c49a:857e5e7908422b6ea5016a3c313f67087fbe2f8b-0' }"

    input:
        tuple val(meta), path(idxml_in), path(mzml)

    output:
        tuple val(meta), path('*.idXML'), emit: idxml
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def prefix = idxml_in.baseName
    def fragment_error = params.fragment_mass_tolerance * 2

    """
    ms2pip_cli.py \\
        --input_idxml $idxml_in \\
        --input_mzml $mzml \\
        --output_idxml ${prefix}_ms2pip.idXML \\
        --num_hits ${params.num_hits} \\
        --model_name ${params.ms2pip_model_name} \\
        --fragment_error $fragment_error \\
        --variable_mods '${params.variable_mods}' \\
        --fixed_mods '${params.fixed_mods}' \\
        --num_cpus ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MS2PIP: \$(conda list | grep "ms2pip" | awk 'NR==2 {print \$2}')
    END_VERSIONS
    """
}

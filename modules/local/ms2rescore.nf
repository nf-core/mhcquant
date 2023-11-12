process MS2RESCORE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ms2rescore=3.0.0b1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ms2rescore:3.0.0b1--pyhdfd78af_1':
        'biocontainers/ms2rescore:3.0.0b1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(idxml), path(mzml), path(fasta)

    output:
    tuple val(meta), path("*ms2rescore.idXML"), emit: rescored_idxml
    path "versions.yml"                       , emit: versions
    // TODO add parsing of the html report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_ms2rescore"
    def ms2pip_model = params.ms2pip_model_name ? "--ms2pip_model ${params.ms2pip_model_name}" : ''
    def ms2_tolerance = 2 * float(params.fragment_mass_tolerance)
    def rescoring_engine = params.rescoring_engine ? "--rescoring_engine ${params.rescoring_engine}"

    """
    ms2rescore.py \\
       --psm_file $idxml \\
       --spectrum_path $mzml \\
       --output_path ${prefix}.idXML \\
       --processes $task.cpus \\
       --feature_generators basic,ms2pip,deeplc \\
       $ms2pip_model \\
       $ms2_tolerance \\
       $rescoring_engine \\
       $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MS²Rescore: \$(echo \$(ms2rescore --version 2>&1) | grep -oP 'MS²Rescore \(v\K[^\)]+' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_ms2rescore"

    """
    touch ${prefix}.idXML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MS²Rescore: \$(echo \$(ms2rescore --version 2>&1) | grep -oP 'MS²Rescore \(v\K[^\)]+' ))
    END_VERSIONS
    """
}

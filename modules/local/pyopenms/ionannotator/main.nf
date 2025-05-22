process PYOPENMS_IONANNOTATOR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyopenms:3.3.0--py313h9b5bd11_0' :
        'biocontainers/pyopenms:3.3.0--py313h9b5bd11_0' }"

    input:
    tuple val(meta), path(mzml), path(fdr_filtered_idxml)

    output:
    tuple val(meta), path("*.tsv")  , emit: tsv
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args  ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def xions  = params.use_x_ions ? "--use_x_ions" : ""
    def zions  = params.use_z_ions ? "--use_z_ions" : ""
    def aions  = params.use_a_ions ? "--use_a_ions" : ""
    def cions  = params.use_c_ions ? "--use_c_ions" : ""

    """
    get_ion_annotations.py \\
        --input $mzml \\
        -idxml $fdr_filtered_idxml \\
        --prefix $meta.id \\
        $args \\
        $xions \\
        $zions \\
        $aions \\
        $cions


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyopenms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_all_peaks.tsv
    touch ${prefix}_matching_ions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyopenms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """
}

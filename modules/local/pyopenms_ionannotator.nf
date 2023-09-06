process PYOPENMS_IONANNOTATOR {
    tag "$sample"
    label 'process_high'

    conda "bioconda::pyopenms=2.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyopenms:2.8.0--py310h3dc0cdb_1' :
        'biocontainers/pyopenms:2.8.0--py310h3dc0cdb_1' }"

    input:
        tuple val(sample), path(mzml), path(fdr_filtered_idxml)

    output:
        tuple val(sample), path("*.tsv"), path("*.tsv"), emit: tsv
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${mzml.baseName}"
        def args             = task.ext.args  ?: ''

        def xions            = params.use_x_ions ? "--use_x_ions" : ""
        def zions            = params.use_z_ions ? "--use_z_ions" : ""
        def aions            = params.use_a_ions ? "--use_a_ions" : ""
        def cions            = params.use_c_ions ? "--use_c_ions" : ""

        """
        get_ion_annotations.py \\
            --input $mzml \\
            -idxml $fdr_filtered_idxml \\
            --prefix $sample \\
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
}

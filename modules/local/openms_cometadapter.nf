process OPENMS_COMETADAPTER {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms-thirdparty:2.6.0--0' :
        'quay.io/biocontainers/openms-thirdparty:2.6.0--0' }"

    input:
        tuple val(meta), path(mzml), path(fasta)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path "versions.yml"             , emit: versions

    script:
        def prefix           = task.ext.prefix ?: "${mzml.baseName}"
        def args             = task.ext.args  ?: ''

        def mods             = params.fixed_mods != " " ? "-fixed_modifications ${params.fixed_mods.tokenize(',').collect { "'${it}'"}.join(" ")}" : ""
        def xions            = params.use_x_ions ? "-use_X_ions true" : ""
        def zions            = params.use_z_ions ? "-use_Z_ions true" : ""
        def aions            = params.use_a_ions ? "-use_A_ions true" : ""
        def cions            = params.use_c_ions ? "-use_C_ions true" : ""
        def nlions           = params.use_NL_ions ? "-use_NL_ions true" : ""
        def remove_precursor = params.remove_precursor_peak ? "-remove_precursor_peak yes" : ""

        """
        CometAdapter -in $mzml \\
            -out ${prefix}.idXML \\
            -database $fasta \\
            -threads $task.cpus \\
            $args \\
            $mods \\
            $xions \\
            $zions \\
            $aions \\
            $cions \\
            $nlions \\
            $remove_precursor

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            openms-thirdparty: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}

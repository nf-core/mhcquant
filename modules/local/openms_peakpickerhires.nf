process OPENMS_PEAKPICKERHIRES {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::openms=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.6.0--h4afb90d_0' :
        'quay.io/biocontainers/openms:2.6.0--h4afb90d_0' }"

    input:
        tuple val(meta), path(mzml)

    output:
        tuple val(meta), path("*.mzML"), emit: mzml
        path "versions.yml"            , emit: versions

    script:
        def prefix           = task.ext.prefix ?: "${mzml.baseName}"

        """
        PeakPickerHiRes -in $mzml \\
            -out ${prefix}.mzML \\
            -algorithm:ms_levels ${params.pick_ms_levels}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}

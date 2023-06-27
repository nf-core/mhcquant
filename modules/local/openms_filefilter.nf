process OPENMS_FILEFILTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::openms=2.9.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.9.1--h135471a_1' :
        'biocontainers/openms:2.9.1--h135471a_1' }"

    input:
        tuple val(meta), path(mzml)

    output:
        tuple val(meta), path("*.mzML"), emit: cleaned_mzml
        path "versions.yml"            , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${mzml.baseName}_cleaned"
        """
        FileFilter -in $mzml \\
            -out ${prefix}.mzML \\
            -peak_options:rm_pc_charge 0 \\
            -threads $task.cpus

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}

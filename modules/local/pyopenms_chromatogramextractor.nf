process PYOPENMS_CHROMATOGRAMEXTRACTOR {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pyopenms=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyopenms:3.1.0--py311h9b8898c_0' :
        'biocontainers/pyopenms:3.1.0--py311h9b8898c_0' }"

    input:
        tuple val(meta), path(mzml)

    output:
        tuple val(meta), path("*.csv")  , emit: csv
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${mzml.baseName}"
        def args             = task.ext.args  ?: ''

        """
        chromatogram_extractor.py \\
            -in $mzml \\
            -out ${prefix}_chrom.csv \\

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pyopenms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${mzml.baseName}"

        """
        touch ${prefix}_chrom.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pyopenms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}

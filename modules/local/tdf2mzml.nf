process TDF2MZML {
    tag "$meta.id"

    container "docker.io/mfreitas/tdf2mzml"

    input:
        tuple val(meta), path(tdf)

    output:
        tuple val(meta), path("*.mzML"), emit: mzml
        path "versions.yml"            , emit: versions

    script:
        def prefix           = task.ext.prefix ?: "${tdf.simpleName}"

        """
        tdf2mzml.py -i $tdf -o ${prefix}.mzML


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version | cut -d ' ' -f2)
            tdf2mzml: \$(echo 0.3.0)
        END_VERSIONS
        """
}

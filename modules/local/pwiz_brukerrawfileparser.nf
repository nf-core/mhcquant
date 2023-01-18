process PWIZ_BRUKERRAWFILEPARSER {
    tag "$meta.id"
    label 'process_long'


    if (params.enable_conda && meta.ext == 'd') {
            { exit 1, "Converting bruker tdf file format to mzml is only supported using docker/singularity. Aborting." }
    }

    container "chambm/pwiz-skyline-i-agree-to-the-vendor-licenses"

    input:
        tuple val(meta), path(tdf)

    output:
        tuple val(meta), path("*.mzML"), emit: mzml
        path "versions.yml"            , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        //def prefix           = task.ext.prefix ?: "${tdf.baseName}"


        """
        wine msconvert \\
            $tdf \\
            -o $PWD \\
            --inten64 \\
            --zlib \\
            --combineIonMobilitySpectra \\
            --filter "scanSumming precursorTol=0.05 scanTimeTol=5 ionMobilityTol=0.1"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            wine: \$(wine --version)
        END_VERSIONS
        """
}

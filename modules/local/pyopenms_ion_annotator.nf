process PYOPENMS_ION_ANNOTATOR {
    tag "$sample"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::pyopenms=2.8.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyopenms:2.8.0--py310h3dc0cdb_1' :
        'quay.io/biocontainers/pyopenms:2.8.0--py310h3dc0cdb_1' }"

    input:
        tuple val(sample), path(mzml), path(fdr_filtered_idxml)

    output:
        tuple val(sample), path("*.tsv"), path("*.tsv"), emit: tsv
        path "versions.yml"             , emit: versions

    script:
        def prefix           = sample
        def precursor_charge = params.prec_charge
        def remove_precursor = params.remove_precursor_peak ? "-remove_precursor_peak" : ""
        def frag_mass_tol    = params.fragment_mass_tolerance
        def xions            = params.use_x_ions ? "-use_x_ions" : ""
        def zions            = params.use_z_ions ? "-use_z_ions" : ""
        def aions            = params.use_a_ions ? "-use_a_ions" : ""
        def cions            = params.use_c_ions ? "-use_c_ions" : ""

        """
        echo $mzml
        get_ion_annotations.py --input $mzml \\
            -idxml $fdr_filtered_idxml \\
            --prefix $sample \\
            --precursor_charge $precursor_charge \\
            --fragment_mass_tolerance $frag_mass_tol \\
            $remove_precursor \\
            $xions \\
            $zions \\
            $aions \\
            $cions \\


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pyopenms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}

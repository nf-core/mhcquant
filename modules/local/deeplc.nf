process DEEPLC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::deeplc=2.2.0 bioconda::pyopenms=2.9.1"
    container 'ghcr.io/jonasscheid/mhcquant:deeplc'

    input:
        tuple val(meta), path(idxml_in)

    output:
        tuple val(meta), path('*.idXML'), emit: idxml
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def prefix = idxml_in.baseName
    def add_abs_rt_error =  params.deeplc_add_abs_rt_error ? "--add_abs_rt_error" : ""
    def add_sqr_rt_error =  params.deeplc_add_sqr_rt_error ? "--add_sqr_rt_error" : ""
    def add_log_rt_error =  params.deeplc_add_log_rt_error ? "--add_log_rt_error" : ""

    """
    deeplc_cli.py \\
        --input $idxml_in \\
        --output ${prefix}_deeplc.idXML \\
        --calibration_mode ${params.deeplc_calibration_mode} \\
        --calibration_bins ${params.deeplc_calibration_bins} \\
        $add_abs_rt_error \\
        $add_sqr_rt_error \\
        $add_log_rt_error

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DeepLC: \$(deeplc --version)
    END_VERSIONS
    """
}

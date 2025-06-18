process SDRF_CONVERT {
    tag "$sdrf"
    label 'process_low'
    
    conda "bioconda::sdrf-pipelines=0.0.20 conda-forge::pandas=1.3.5 python=3.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sdrf-pipelines:0.0.20--pyhdfd78af_0' :
        'biocontainers/sdrf-pipelines:0.0.20--pyhdfd78af_0' }"
        
    input:
    path sdrf
    
    output:
    path "samplesheet.tsv", emit: samplesheet
    path "versions.yml", emit: versions
    
    script:
    """
    python ${projectDir}/bin/sdrf_to_mhcquant.py \\
        -s $sdrf \\
        -o samplesheet.tsv
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sdrf-pipelines: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('sdrf-pipelines').version)")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
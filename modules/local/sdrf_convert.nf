process SDRF_CONVERT {
    tag "$sdrf"
    label 'process_low'
    
    conda "bioconda::sdrf-pipelines=0.0.20"
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
    parse_sdrf.py convert-mhcquant \\
        -s $sdrf \\
        -o samplesheet.tsv
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sdrf-pipelines: \$(parse_sdrf.py --version | sed 's/parse_sdrf.py version //g')
    END_VERSIONS
    """
}
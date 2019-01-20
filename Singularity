From:nfcore/base
Bootstrap:docker

%labels
    DESCRIPTION Singularity image containing all requirements for the nf-core/mhcquant pipeline
    VERSION 1.2.0

%environment
    PATH=/opt/conda/envs/nf-core-mhcquant-1.2.0/bin:$PATH
    PATH=/opt/conda/envs/nf-core-mhcquant-percolator-1.0dev/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda env create -f /environment-percolator.yml
    /opt/conda/bin/conda clean -a

FROM nfcore/base
MAINTAINER Leon Bichmann <leon.bichmann@informatik.uni-tuebingen.de>
LABEL description="Docker image containing all requirements for nf-core/mhcquant pipeline"

COPY environment.yml /
COPY environment-percolator.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env create -f /environment-percolator.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-mhcquant-1.0.0/bin:$PATH
ENV PATH /opt/conda/envs/nf-core-mhcquant-percolator-1.0dev/bin:$PATH

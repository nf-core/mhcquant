FROM nfcore/base
LABEL authors="Leon Bichmann" \
      description="Docker image containing all requirements for nf-core/mhcquant pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env create -f /environment-percolator.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-mhcquant-1.2.7dev/bin:$PATH
ENV PATH /opt/conda/envs/nf-core-mhcquant-percolator-1.2.6/bin:$PATH

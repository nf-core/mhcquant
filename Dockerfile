FROM nfcore/base
MAINTAINER Leon Bichmann <leon.bichmann@informatik.uni-tuebingen.de>
LABEL description="Docker image containing all requirements for nf-core/openmspeptidequant pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-openmspeptidequant-1.0dev/bin:$PATH


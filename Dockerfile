FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/openmspeptidequant pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN mkdir /opt/thirdparty && wget https://github.com/OpenMS/THIRDPARTY/raw/master/Linux/64bit/Comet/comet.exe && chmod ugo+x comet.exe && mv comet.exe /opt/thirdparty
ENV PATH /opt/conda/envs/nf-core-openmspeptidequant-1.0dev/bin:$PATH
ENV PATH /opt/conda/envs/nf-core-openmspeptidequant-1.0dev/opt/thirdparty:$PATH


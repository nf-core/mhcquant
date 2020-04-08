FROM nfcore/base:1.9
LABEL authors="Leon Bichmann, Lukas Heumos, Alexander Peltzer" \
      description="Docker image containing all requirements for nf-core/mhcquant pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-mhcquant-1.5dev > nf-core-mhcquant-1.5dev.yml

# Add conda installation dir and percolator to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-mhcquant-1.5dev/bin:$PATH

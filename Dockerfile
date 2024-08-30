#
# MM30 GEO Data Docker Image Construction
# V. Keith Hughitt
#
FROM rocker/r-ver:4
MAINTAINER keith.hughitt@nih.gov

# conda/snakemake setup
ENV PATH /opt/conda/bin:${PATH}
ENV LANG C.UTF-8
ENV SHELL /bin/bash

RUN apt-get update

RUN apt-get install -y wget curl bzip2 ca-certificates gnupg2 squashfs-tools git

RUN /bin/bash -c "curl -L https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh > miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh"

RUN /bin/bash -c "conda install -y -c conda-forge mamba && \
    mamba create -q -y -c conda-forge -c bioconda -n snakemake && \
    source activate snakemake && \
    mamba install -q -y -c conda-forge -c bioconda \
    snakemake-minimal frictionless r-annotables r-arrow r-r.utils r-tidyverse bioconductor-geoquery bioconductor-biomart"

RUN echo "source activate snakemake" > ~/.bashrc
ENV PATH /opt/conda/envs/snakemake/bin:${PATH}

# copy code over
WORKDIR /geo

COPY R/ R/
COPY python/ python/
COPY Snakefile Snakefile
COPY metadata/ metadata/
COPY supp/ supp/

# launch bash when container is started
ENTRYPOINT ["bash"]

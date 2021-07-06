#
# MM25 GEO Data Docker Image Construction
# V. Keith Hughitt
#
# - built on top of the rocker/4.1.0 image
# - renv is used to retrieve specific package versions
#
FROM rocker/r-ver:4.1.0
MAINTAINER keith.hughitt@nih.gov

# install libcurl
RUN apt update
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libpng-dev libxml2-dev

# install renv
ENV RENV_VERSION 0.13.2
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# copy code over
WORKDIR /geo
COPY annot/ annot/
COPY extra/ extra/
COPY R/ R/
COPY supp/ supp/
COPY renv.lock renv.lock

# restore renv snapshot
RUN R -e 'renv::restore()'

# launch bash when container is started
ENTRYPOINT ["bash"]

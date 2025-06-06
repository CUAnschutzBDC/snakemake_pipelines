# Get base image
FROM bioconductor/bioconductor_docker:RELEASE_3_21

# environment
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV RENV_VERSION=1.0.7

# Install R packages with renv
RUN R -e "install.packages('renv', version = '${RENV_VERSION}')"


WORKDIR /home/rstudio
COPY renv.lock renv.lock
COPY R_dependencies.R R_dependencies.R
RUN R -e "renv::restore()"
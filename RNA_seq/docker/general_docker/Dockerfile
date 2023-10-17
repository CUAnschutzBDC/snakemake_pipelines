# Get base image
FROM bioconductor/bioconductor_docker:RELEASE_3_16

# environment
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH $PATH:/opt/conda/bin
ENV RENV_VERSION 0.16.0

# get conda installed
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean --all -f -y && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    apt-get update

COPY environment.yaml .
RUN /opt/conda/bin/conda config --add channels defaults && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda config --add channels agbiome && \
    /opt/conda/bin/conda install -n base conda-libmamba-solver && \
    conda config --set solver libmamba && \
    conda env update -n base -f environment.yaml && \
    conda clean -afy

# Install R packages with renv
#RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
#RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"


#WORKDIR /home/rstudio
#COPY renv.lock renv.lock
#COPY R_dependencies.R R_dependencies.R
#RUN R -e "renv::restore()"
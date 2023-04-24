ARG UBUNTU_VER=18.04
#ARG CONDA_VER=latest
ARG CONDA_VER=4.7.12.1
ARG OS_TYPE=x86_64
#ARG PY_VER=3.9.0
ARG SM_VER=5.26.1

FROM ubuntu:${UBUNTU_VER}

ADD . /home/qcgenes
WORKDIR /home/qcgenes

# System packages
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -yq curl wget jq vim git libtbb2 zip tzdata

# Use the above args during building https://docs.docker.com/engine/reference/builder/#understand-how-arg-and-from-interact
ARG CONDA_VER
ARG OS_TYPE
# Install miniconda to /miniconda
RUN curl -LO "http://repo.continuum.io/miniconda/Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh"
RUN bash Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh -p /miniconda -b
RUN rm Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda init bash
RUN conda update -y conda

#ARG PY_VER
ARG SM_VER
#RUN conda install -c anaconda -y python=${PY_VER}
RUN conda install -c conda-forge -y mamba
RUN mamba install -c conda-forge -c bioconda -y snakemake=${SM_VER}
RUN snakemake --use-conda --conda-frontend mamba --conda-create-envs-only --cores 1


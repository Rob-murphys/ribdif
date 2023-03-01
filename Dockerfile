# Dockerfile, Image, Container

FROM continuumio/miniconda3

MAINTAINER RobMur

RUN conda create -n ribdif python=3.11 -y \
    && echo "source activate ribdif" > ~/.bashrc \
    && conda install -c bioconda vsearch fasttree mafft barrnap pyani -y \
    && pip install --no-cache-dir pandas==1.5.2 biopython==1.80 numpy==1.24.1 matplotlib==3.6.2 seaborn==0.12.2 scipy==1.10.0 ncbi_genome_download==0.3.1 networkx==3.0 fastcluster==1.2.6 chardet==5.1.0

ENV PATH /opt/conda/envs/env/bin:$PATH

WORKDIR /ribdif

COPY ribdif/ ./

COPY docs/default.primers ./default.primers

ENTRYPOINT [ "python", "__main__.py" ]



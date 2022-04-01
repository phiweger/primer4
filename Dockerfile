FROM continuumio/miniconda3:4.10.3
ENV VERSION 3.8
ENV TOOL python
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update --allow-releaseinfo-change && apt install -y procps wget gzip git pigz bc && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels default

# And some more dependencies
RUN conda install -y python=3.8 pandas=1.3.5 numpy=1.22.0 click=8.0.3 pip=21.3.1 pytest=6.2.5 tqdm=4.62.3

RUN conda install -y -c bioconda hgvs=1.5.1 blast=2.12.0 pyfaidx=0.6.3.1 gffutils=0.10.1 pysam=0.17.0 primer3-py=0.6.1

RUN pip install cdot==0.2.2
RUN git clone https://github.com/phiweger/primer4 && cd primer4 && pip install -e .
# hgvs alternative https://github.com/counsyl/hgvs

RUN conda clean -y -a

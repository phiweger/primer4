<p align="center">
    <img src="img/logo.png" alt="Logo" width="700">
</p>


## README

_Please note: This repo is in full development, please handle with care._

In human genetics, an important task is validating mutations that have been found using exome or whole genome sequencing with a PCR spanning the abnormal region. This task is not difficult, but surprisingly laborious when done manually.

This workflow automates the design of PCR and qPCR primers. As input, it takes either:

- an HGVS-formatted mutation, e.g. `NM_000546.6:c.215C>G` (for PCR)
- a gene name and exon number, e.g. `NM_007294.4::14` (for qPCR)

The workflow will:

- validate [HGVS](https://varnomen.hgvs.org/bg-material/simple/) syntax
- for mutations in coding regions try to span the exon
- constrain primers to desired parameters
- construct a second backup pair of primers that does not overlap the first pair
- avoid placing primers across known common variants
- validate all primers to not produce unspecific amplicons


## Setup

```bash
NAME=primer4
conda create -y -n $NAME -c bioconda -c conda-forge mamba python=3.8
conda activate $NAME

mamba install -y -c bioconda -c conda-forge \
    pandas \
    numpy \
    click \
    pip \
    pytest \
    streamlit=1.8.1 \
    hgvs=1.5.1 \
    blast=2.12.0 \
    pyfaidx=0.6.3.1 \
    gffutils=0.10.1 \
    pysam=0.17.0 \
    primer3-py=0.6.1 \
    pygenometracks=3.7

pip install cdot==0.2.2
# streamlit-aggrid==0.3.3

# conda list -n $NAME | grep -v '#' > requirements.txt

git clone https://github.com/phiweger/primer4 && cd primer4 && pip install -e .
```


## Test

```bash
pytest --config config.json
```


## Run

```bash
streamlit run primer4/stream.py -- --config config.json
# Args after "--":
# https://discuss.streamlit.io/t/command-line-arguments/386/4
```

You should see something like this (modify `tracks.empty.ini` to adjust the plot; additional options are documented in `tracks.annotation.ini`):

<p align="center">
    <img src="img/interface.png" alt="Interface" width="700">
</p>

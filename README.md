## README

In human genetics, an important task is validating mutations that have been found using exome or whole genome sequencing with a PCR spanning the abnormal region. This task is not difficult, but surprisingly laborious when done manually.

This workflow automates the design of PCR and qPCR primers. As input, it takes either:

- an HGVS-formatted mutation, e.g. `NM_000546.6:c.215C>G` (for PCR)
- a gene name and exon number, e.g. `NM_007294.4:14` (for qPCR)

The workflow will:

- validate [HGVS](https://varnomen.hgvs.org/bg-material/simple/) syntax
- for mutations in coding regions try to span the exon
- constrain primers to desired parameters
- construct a second backup pair of primers that does not overlap the first pair
- avoid placing primers across known common variants
- validate all primers to not produce unspecific amplicons


## Setup

```bash
conda install -y -n primer -c bioconda nextflow && conda activate primer
git clone https://github.com/phiweger/primer && cd primer
mkdir data
# ^ add data in there, see below

# And some more dependencies
conda install -y python=3.8 pandas=1.3.5 numpy=1.22.0 click=8.0.3 pip=21.3.1 pytest=6.2.5

conda install -y -c bioconda hgvs=1.5.1 blast=2.12.0 pyfaidx=0.6.3.1 gffutils=0.10.1 pysam=0.17.0 primer3-py=0.6.1

pip install cdot=0.2.2
git clone github.com/phiweger/primer4 && cd primer4 && pip install -e .
# hgvs alternative https://github.com/counsyl/hgvs
```

Then set the absolute (!) path to the data in `nextflow.config`. If you use other reference genomes/ annotations/ variant databases, you also have to adjust the filenames in `config.json`.



## Run

```bash
primer4 -i order.csv -d /path/to/primer4/data -p /path/to/primer4/config.json

# DEPRECATED
cat input.csv
# name,method,variant
# ABCA4_v1,PCR,NM_000350.3:c.4234C>T
# ABCA4_v2,PCR,NM_000350.3:c.4773+3A>G

nextflow run workflow/main.nf --input input.csv --results designs
```


## Run with Docker

```bash
docker build -t foobar .
# This part is important! .../primer4/data:/primer4/data
docker run -v $PWD:/workdir -v .../primer4/data:/primer4/data -it -t foobar /bin/bash
primer4 -i /workdir/order.csv -d /primer4/data/ -p /workdir/config.json
```



## Data provenance

You need to get three things, which need to correspond to one another:

1. reference genome
2. SNPs
3. annotation

Below, we'll use the files corresponding to the human genome `hg19` (v13, [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/), last access 2022-02-03).


```bash
# 1. Reference genome and index for random access:
# - GRCh37_latest_genomic.fna.gz
# - GRCh37_latest_genomic.fna.gz.fai

wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
gunzip GRCh37_latest_genomic.fna.gz
samtools faidx GRCh37_latest_genomic.fna

# Blast index
makeblastdb -in GRCh37_latest_genomic.fna -dbtype nucl
# Adding sequences from FASTA; added 297 sequences in 39.731 seconds.
```


```bash
# 2. Variants and index:
# - GRCh37_latest_dbSNP_all.vcf.gz
# - GRCh37_latest_dbSNP_all.vcf.gz.tbi

wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_dbSNP_all.vcf.gz
gunip GRCh37_latest_dbSNP_all.vcf.gz
bgzip -c GRCh37_latest_dbSNP_all.vcf > GRCh37_latest_dbSNP_all.vcf.gz
tabix -p vcf GRCh37_latest_dbSNP_all.vcf.gz

# To split and concat
# split -b 4000M -d GRCh37_latest_dbSNP_all.vcf.gz vcf
# cat vcf* > GRCh37_latest_dbSNP_all.vcf.gz
osf -p bc7sr clone
cat vcf* > GRCh37_latest_dbSNP_all.vcf.gz

# https://genetools.org/SNPCheck/snpcheck.htm
DB=ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf
gunzip ${DB}.gz
bgzip -c $DB > ${DB}.gz
tabix -p vcf ${DB}.gz

DB=ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf
tar -xzvf ${DB}.tar.gz
cat *.vcf > $DB
bgzip -c $DB > ${DB}.gz
tabix -p vcf ${DB}.gz
```


```bash
# 3. Annotation database:
# - hg19-p13_annotation.db
```

Get it from [OSF](https://osf.io) (project ID [7csav](https://osf.io/7csav/)); here is how we created it:


```bash
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz
gunzip GRCh37_latest_genomic.gff.gz
# switch to Python
```

```python
import gffutils
# http://daler.github.io/gffutils/
# https://daler.github.io/gffutils/database-import.html

fp = 'GRCh37_latest_genomic.gff'
db = gffutils.create_db(
    fp, 
    dbfn='hg19-p13_annotation.db',
    force=True,
    keep_order=True,
    merge_strategy='merge',
    sort_attribute_values=True)
# Takes a couple of hours on a regular laptop, then use like
db = gffutils.FeatureDB('hg19-p13_annotation.db', keep_order=True)
db['exon-NR_024540.1-3']
```


## Design choices

For future reference, below are documented some peculiarities:

- `gffutils.FeatureDB` will not accept symlinks; this means we cannot pass the database in a `nextflow` channel. 
- There are two files with parameters, one for `nextflow` (`nextflow.config`) and one to guide the primer design (`config.json`). There are more elegant solutions, but meh.

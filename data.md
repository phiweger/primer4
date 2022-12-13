## Data

Get it from [OSF](https://osf.io) (project ID [7csav](https://osf.io/7csav/)); below you find how we created them and how you find the data. Before running `primer4`, you need to set the paths in the `config.json` accordingly (note that versions need to correspond to one another):

1. (human) reference genome
2. transcript annotation
3. variants
4. chromosome name glossary
5. transcript coordinate mappings
6. transcript sequence data < optional

Yeah, you only wanted to design some primers ...

---

- GRCh37, deprecated, last accessed 2022-02-03
- GRCh38, last accessed 2022-11-21, https://www.ncbi.nlm.nih.gov/genome/guide/human/

```bash
mkcd 4primer4
conda activate primer4

# Human reference genome
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz
samtools faidx GRCh38_latest_genomic.fna
makeblastdb -in GRCh38_latest_genomic.fna -dbtype nucl

# Transcripts
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
gunzip GRCh38_latest_genomic.gff.gz
```

---

To query the annotation during primer search, we need to index and organize them
using `gffutils`. This takes forever, but is only done once per genome release.


```python
import gffutils
# http://daler.github.io/gffutils/
# https://daler.github.io/gffutils/database-import.html
# https://daler.github.io/gffutils/autodocs/gffutils.create.create_db.html
# https://github.com/daler/gffutils/issues/20
# https://www.biostars.org/p/152517/#152604
# https://github.com/daler/gffutils/issues/48
# https://daler.github.io/gffutils/database-import.html
# https://daler.github.io/gffutils/autodocs/gffutils.create.create_db.html

db = gffutils.create_db(
    fp,
    dbfn='hg38-p14_annotation.db',
    force=True,
    keep_order=True,
    merge_strategy='merge',
    sort_attribute_values=True)
# 2 days, 3.5 gb
```

---

We also need to create a glossary of chromosome names, because
different variant databases use different names (of course).

```bash
pip install screed
python ../scripts/chrom_names.py --genome GRCh37_latest_genomic.fna --out chrom_names_hg38.csv
```

---

To query SNVs on primers, we need a couple of databases, which are selected after the [`snpcheck`](https://genetools.org/SNPCheck/docs.htm) webservice.

```bash
# https://genetools.org/SNPCheck/snpcheck.htm

# dbSNP (same site as reference genome and annotation above ^)
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_dbSNP_all.vcf.gz
gunzip GRCh38_latest_dbSNP_all.vcf.gz
bgzip -c GRCh38_latest_dbSNP_all.vcf > GRCh38_latest_dbSNP_all.vcf.gz
tabix -p vcf GRCh38_latest_dbSNP_all.vcf.gz
# To split and concat
# split -b 4000M -d GRCh37_latest_dbSNP_all.vcf.gz vcf
# cat vcf* > GRCh37_latest_dbSNP_all.vcf.gz
osf -p bc7sr clone
cat vcf* > GRCh37_latest_dbSNP_all.vcf.gz

# 1000 genomes project:
# https://www.internationalgenome.org/data
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz
DB=ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf
gunzip ${DB}.gz
bgzip -c $DB > ${DB}.gz
tabix -p vcf ${DB}.gz

# https://evs.gs.washington.edu/EVS/ and go to "Downloads" tab
wget http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz
DB=ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf
tar -xzvf ${DB}.tar.gz
cat *.vcf > $DB
bgzip -c $DB > ${DB}.gz
tabix -p vcf ${DB}.gz
```

Note on the latter:

> The bulk files of the ESP 6500 exome data below are still primarily in GRCh37 (or HG19), the GRCh38 lifted-over positions are added in an extra column in the text file, or in an extra attribute in the INFO field in the VCF file.

---

To parse the variants correctly and annotate their position on the transcript, we use `cdot`.

```bash
# Transcript coordinate mappings:
# - cdot-0.2.1.refseq.grch37_grch38.json.gz
# cdot replaces UTA lookups from hgvs library
# https://github.com/SACGF/cdot

# Get it from https://osf.io/7csav/
# Alternative (check version):
# wget https://cdot.cc/download/cdot-0.2.8.refseq.grch37_grch38.json.gz
```

The following data is optional. If `primer4` cannot find it locally, it will call an API (ie you need internet access when using the tool).

```bash
# Transcript sequence data < optional
# For local execution:
# - https://github.com/biocommons/hgvs/issues/634
# - https://hgvs.readthedocs.io/en/stable/installation.html#installing-seqrepo-optional
# There is even a docker container
# https://github.com/biocommons/biocommons.seqrepo/blob/main/docs/docker.rst
pip install biocommons.seqrepo==0.6.5
sudo seqrepo -r seqrepo pull
# Then management is taken care of by primer4, ie this is not necessary:
# export HGVS_SEQREPO_DIR=${PWD}/seqrepo/2021-01-29
```

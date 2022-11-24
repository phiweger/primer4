## Data

GRCh38, last accessed 2022-11-21, https://www.ncbi.nlm.nih.gov/genome/guide/human/


```bash
mkcd 4primer4
conda activate primer4

# Human reference genome
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz
samtools faidx GRCh38_latest_genomic.fna
makeblastdb -in GRCh38_latest_genomic.fna -dbtype nucl

# Variants
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_dbSNP_all.vcf.gz
gunzip GRCh38_latest_dbSNP_all.vcf.gz
bgzip -c GRCh38_latest_dbSNP_all.vcf > GRCh38_latest_dbSNP_all.vcf.gz
tabix -p vcf GRCh38_latest_dbSNP_all.vcf.gz

# Transcripts
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
gunzip GRCh38_latest_genomic.gff.gz
```


```python
import gffutils
# http://daler.github.io/gffutils/
# https://daler.github.io/gffutils/database-import.html
# https://daler.github.io/gffutils/autodocs/gffutils.create.create_db.html
# https://github.com/daler/gffutils/issues/20
# https://www.biostars.org/p/152517/#152604
# https://github.com/daler/gffutils/issues/48


fp = 'GRCh38_latest_genomic.gff'
# start 20:57
db = gffutils.create_db(
    fp, 
    dbfn='hg38-p14_annotation.db',
    force=True,
    keep_order=True,
    merge_strategy='merge',
    sort_attribute_values=True)


db = gffutils.create_db(
    fp, 
    dbfn='hg38-p14_annotation.db',
    force=True,
    keep_order=True,
    merge_strategy='merge',
    sort_attribute_values=True,
    verbose=True)

db = gffutils.create_db(fp, dbfn='hg38-p14_annotation.db', verbose=True)

import gffutils
# http://daler.github.io/gffutils/
# https://daler.github.io/gffutils/database-import.html
# https://daler.github.io/gffutils/autodocs/gffutils.create.create_db.html


fp = 'GRCh38_latest_genomic.gff'
```




In [3]: db = gffutils.create_db(
   ...:     fp,
   ...:     dbfn='hg38-p14_annotation.db',
   ...:     force=True,
   ...:     keep_order=True,
   ...:     merge_strategy='merge',
   ...:     sort_attribute_values=True)

2 days, 3.5 gb


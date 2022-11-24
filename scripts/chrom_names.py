import screed
import re


fp = 'GRCh37_latest_genomic.fna'
pattern = r'.*chromosome ([XY0-9]{1,})[ ,].*'
'''
NC_000021.8 Homo sapiens chromosome 21, GRCh37.p13 Primary Assembly
NT_113950.2 Homo sapiens chromosome 21 unlocalized genomic contig, GRCh37.p13 Primary Assembly
'''

d = {}
with screed.open(fp) as file:
    for line in file:
        if (not 'unplaced' in line.name) and (not 'mitochondrion' in line.name):
            x = re.match(pattern, line.name).group(1)
            d[line.name.split(' ')[0]] = x


with open('chrom_names_hg38.csv', 'w+') as out:
    for k, v in d.items():
        out.write(f'{k},{v}\n')


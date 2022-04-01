#! /usr/bin/env python3


from collections import defaultdict
import json
import pdb
from pathlib import Path

import click
# https://click.palletsprojects.com/en/8.0.x/arguments/
import gffutils
# http://daler.github.io/gffutils/
# import primer3
from pyfaidx import Fasta
from pysam import VariantFile
# https://pysam.readthedocs.io/en/latest/usage.html#working-with-vcf-bcf-formatted-files
from screed import rc
# from tqdm import tqdm

from primer.utils import (
    infer_coordinates,
    variant_context,
    common_variants,
    design_primers,
    parse_design,
    log)


@click.command()
@click.option('-n', '--pname', default='primer', help='Name for primer pair')
@click.option('-p', '--prefix', default='primer', help='Outfile prefix')
@click.option('-q', '--query', required=True, help='Variant in HGVS syntax (e.g. "NM_000546.6:c.215C>G")')
@click.option('-m', '--method', default='PCR', show_default=True, help='Method primers will be used with', type=click.Choice(['PCR', 'qPCR'], case_sensitive=False))
@click.option('-c', '--config', required=True, help='Path to settings', type=click.Path(exists=True, resolve_path=True))
@click.option('-d', '--data', required=True, help='Path to data', type=click.Path(exists=True, resolve_path=True))
def get_me_some_primers(pname, prefix, query, method, config, data):
    '''
    Usage:

    pythonprimer.py -n zebra -p primers -v 'NM_000546.5:xyz.215C>G' -m 'PCR' -c config.json

    Examples:
    
    NM_000546.6:c.215C>G
    NM_138459.3:c.26G>A
    NM_000492.4:c.1679+34G>T
    NM_000492.4:c.2989-313A>T
    '''
    assert not ' ' in pname, 'No space allowed in primer names, sorry'
    
    
    # 0. Housekeeping
    # PCR params according to "AA Primerdesign", roxtra, last access 2022-01-13
    with open(config, 'r') as file:
        params = json.load(file)

    data = Path(data)

    fp_annotation = str(data / params['annotation'])
    fp_genome = str(data / params['reference'])
    fp_variants = str(data / params['variation'])

    db = gffutils.FeatureDB(fp_annotation, keep_order=True)
    genome = Fasta(fp_genome)
    variants = VariantFile(fp_variants)
    

    # 1. Validate variant conforms to HGVS nomenclature, then
    # 2. Translate to genomic coordinates
    name, chromosome, g_pos = infer_coordinates(query, db)
    # name is tx name, eg NM_000546.6 when query is NM_000546.6:c.215C>G
    

    # 3. Our PCR ideally spans the exon that contains the variant, but which exon?
    template, boundaries, constraints = variant_context(
        name, genome, chromosome, g_pos, db, params)
    
    
    # 4. Mask "common" variants in target region
    '''
    zcat < GRCh37_latest_dbSNP_all.vcf.gz | grep 'COMMON' | head -n1 
    ##INFO=<ID=COMMON,Number=0,Type=Flag,Description="RS is a common SNP.  A common SNP is one that has at least one 1000Genomes population with a minor allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.">
    '''
    masked = common_variants(template, chromosome, boundaries, variants)
    
    
    # 5. Design primers
    design = design_primers(method, params, masked, constraints)
    # pdb.set_trace()
    primers = parse_design(
        design, params['n_return'], f'{pname}-round1', g_pos, chromosome, masked)


    # Now we want to mask all positions against which the best primer has been
    # designed to get another "unrelated" pair. Note: It might still be that
    # the first primer has to be discarded due to non-unique binding.
    first = primers[f'{pname}-round1::c0']

    p1 = range(first['fwd']['start'], first['fwd']['end'])  # not + 1
    p2 = range(first['rev']['start'], first['rev']['end'])  # not + 1

    from itertools import chain
    masked2 = ''.join(
        ['N' if ix in chain(p1, p2) else i for ix, i in enumerate(masked)])

    print(masked)
    print(masked2)
    print(first['fwd']['sequence'])
    print(first['rev']['sequence'])
    
    design2 = design_primers(method, params, masked2, constraints)
    primers2 = parse_design(
        design2, params['n_return'], f'{pname}-round2', g_pos, chromosome, masked2)

    print('\n' + log(f'Found {len(primers)} primer(s)'))
    

    with open(f'{prefix}.template.fna', 'w+') as out:
        out.write(f'>template\n{template}\n')


    with open(f'{prefix}.json', 'w+') as out:
        json.dump({**primers, **primers2}, out, indent=4, sort_keys=True)
    
    
    with open(f'{prefix}.tsv', 'w+') as out:
        for k, v in primers.items():
            out.write(f"{k}\t{v['fwd']['sequence']}\t{v['rev']['sequence']}\n")
        for k, v in primers2.items():
            out.write(f"{k}\t{v['fwd']['sequence']}\t{v['rev']['sequence']}\n")


    # with open(f'{prefix}.fna', 'w+') as out:
    #     for k, v in primers.items():
    #         out.write(f">{k}-fwd\n{v['fwd']['sequence']}\n")
    #         out.write(f">{k}-rev\n{v['fwd']['sequence']}\n")


    print(log('Done. May the thermodynamic force be with you.'))


if __name__ == '__main__':
    get_me_some_primers()


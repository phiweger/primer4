import json
from pathlib import Path
import pdb

from cdot.hgvs.dataproviders import JSONDataProvider
import gffutils
from pyfaidx import Fasta

import click

from primer4.models import Variant, ExonDelta, SingleExon, ExonSpread, Template
from primer4.design import design_primers, check_for_multiple_amplicons
from primer4.utils import mask_sequence, reconstruct_mrna


def gimme_some_primers(method, code, fp_genome, genome, hdp, db, vardbs, params):
    '''
    fp_data = '/Users/phi/Dropbox/repos/primer4/data'
    fp_config = '/Users/phi/Dropbox/repos/primer4/config.json'

    gimme_some_primers('sanger', 'NM_000546.6:c.215C>G', fp_data, fp_config)
    gimme_some_primers('qpcr', ('NM_001145408.2', 6), ...)
    gimme_some_primers('mrna', ('NM_000546.6', 6, 7), ...)
    '''
    if method == 'sanger':
        v = Variant(code[0], hdp, db)
    elif method == 'qpcr':
        v = SingleExon(code[0], int(code[1]))
    elif method == 'mrna':
        v = ExonSpread(code[0], int(code[1]), int(code[2]))
    else:
        raise ValueError('Method is not implemented, exit.')

    tmp = Template(v, db)
    tmp.load_variation_(vardbs)
    
    if method == 'sanger':
        masked = mask_sequence(tmp.get_sequence(genome), tmp.mask)
        constraints = tmp.apply(method, db, params)
        primers = [p for p in next(design_primers(masked, constraints, params, []))]

    elif method == 'qpcr':
        masked = mask_sequence(tmp.get_sequence(genome), tmp.mask)
        
        primers = []
        for constraints in tmp.apply('qpcr', db, params):
            # print(constraints)
            x = [p for p in next(design_primers(masked, constraints, params, []))]
            primers.extend(x)

    elif method == 'mrna':
        tmp.mrna = reconstruct_mrna(tmp.feat, db, genome, vardbs)
        constraints = tmp.apply('mrna', db, params)
        masked = tmp.mrna[0]
        primers = [p for p in next(design_primers(masked, constraints, params, []))]
    
    else:
        raise ValueError('Method is not implemented, exit.')
    
    results = check_for_multiple_amplicons(primers, fp_genome)
    return results


# https://pymbook.readthedocs.io/en/latest/click.html
# https://click.palletsprojects.com/en/8.0.x/options/
@click.command()
@click.option('-i', '--order', type=click.Path(exists=True))
@click.option('-d', '--fp-data', type=click.Path(exists=True), default='primer4/data')
@click.option('-p', '--fp-config', type=click.Path(exists=True), default='config.json')
def main(order, fp_data, fp_config):
    
    print('Housekeeping ...')
    fp_genome = f'{fp_data}/GRCh37_latest_genomic.fna'
    fp_coords = f'{fp_data}/cdot-0.2.1.refseq.grch37_grch38.json.gz'
    fp_annotation = f'{fp_data}/hg19-p13_annotation_bak.db'
    # fp_annotation = f'{fp_data}/hg19-p13_annotation.db'
    
    fp_snvs_1 = f'{fp_data}/GRCh37_latest_dbSNP_all.vcf.gz'
    fp_snvs_2 = f'{fp_data}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz'
    fp_snvs_3 = f'{fp_data}/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.gz'

    genome = Fasta(fp_genome)
    hdp = JSONDataProvider([fp_coords])
    db = gffutils.FeatureDB(fp_annotation, keep_order=True)
    
    with open(fp_config, 'r') as file:
        params = json.load(file)
    
    vardbs = {
        'dbSNP': fp_snvs_1,
        '1000Genomes': fp_snvs_2,
        'ESP': fp_snvs_3
        }

    l = []
    with open(order, 'r') as file:
        for line in file:
            method, code = line.strip().split(',')
            code = code.split('::')

            method = method.lower()

            print('\n\n', method, code)
            primers = gimme_some_primers(method, code, fp_genome, genome, hdp, db, vardbs, params)
            
            for pair in primers:
                print(pair)

            # Path('results').mkdir()
            for pair in primers:
                # if code[0] == 'NM_006087.4:c.745G>A' and pair.penalty < 1:
                    # pdb.set_trace()

                with open(f'results/{pair.name}.json', 'w+') as out: 
                    json.dump(pair.data, out, indent=4, sort_keys=True)
                # print(pair.name)
                # print(pair.data)

                l.append(f'{line.strip()},{pair.name},{pair.penalty}\n')

    with open('results.csv', 'w+') as out:
        for i in l:
            out.write(i)

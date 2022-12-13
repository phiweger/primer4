'''
Remove duplicate sequences from reference genome which might cause Blast to erroneously identify multimappings during primer validation.
'''
import argparse

import screed


parser = argparse.ArgumentParser()
parser.add_argument(
    '--genome', default='GRCh37_latest_genomic.fna', required=True,
    help='Reference genome from NCBI')
parser.add_argument(
    '--out', default='4blast.fna', help='Name outfile')  
args = parser.parse_args()


exclude = set([
    'unlocalized',
    'genomic patch of type FIX',
    'alternate locus group',
    ])

with screed.open(args.genome) as file, open(args.out, 'w+') as out:
    for line in file:
        if any([i in line.name for i in exclude]):
            continue
        else:
            out.write(f'>{line.name}\n{line.sequence}\n')

#! /usr/bin/env python3


from collections import defaultdict
import csv
import json
import pdb

import click
# https://click.palletsprojects.com/en/8.0.x/arguments/
import screed
import pandas as pd


def m13(method='PCR', orient=None):
    if method == 'PCR':
        return 'ohne M13'
    elif method == 'qPCR':
        assert orient, f'Need an orientation with {method}'
        if orient == 'fwd':
            return 'M13_F'
        elif orient == 'rev':
            return 'M13_R'
        else:
            raise ValueError(f'What orientation is "{orient}"?')
    else:
        raise ValueError('Method not implemented')


def smallest_penalty(primers):
    penalty = 1e10
    best = ''
    
    for name, spec in primers.items():
        if spec['penalty'] < penalty:
            penalty = spec['penalty']
            best = name

    return best


@click.command()
@click.option('-c', '--candidates', required=True, help='Primer suggestions')
@click.option('-a', '--annealing', required=True, help='in silico PCR result')
@click.option('-m', '--method', default='PCR', show_default=True, help='Method primers will be used with', type=click.Choice(['PCR', 'qPCR'], case_sensitive=False))
@click.option('-o', '--outfile', required=True, help='Results file')
def validate_primers(candidates, annealing, method, outfile):


    with open(candidates, 'r') as file:
        cand = json.load(file)
    
    
    regions = defaultdict(list)
    with screed.open(annealing) as file:
        for i in file:
            region, name, *rest = i.name.split(' ')
            # Splitting strings on space is a sure footgun in bioinformatics, so
            # just break if there is something unexpected here.
            assert len(rest) + 2 == 5, 'There is a problem with naming ...'
            regions[name].append(region) 
    
    
    l = []
    header = 'version orientation M13 name from to sequence annealing amplicon snv pseudo comment lid box date'.split(' ')

    for name, coords in regions.items():
       # Accept only ...
       if (len(coords) == 1) and (name in cand.keys()):  
       # Predicted to generate a single amplicon
           try:
               spec = cand[name]
           except KeyError:
               continue

           chromosome, startend = coords[0].split(':')
           start, end = startend.split('+')
             
           
           for orient in ['fwd', 'rev']:
               l.append([
                   '', 
                   'F' if orient == 'fwd' else 'R',
                   m13(method, orient),
                   name,
                   start,
                   end,
                   spec[orient]['sequence'],
                   spec[orient]['Tm'],
                   spec['insert'],
                   'ok',
                   'ok',
                   '',
                   '',
                   '',
                   '',
               ])


    # The followign primers do not have secondary alignments
    df = pd.DataFrame.from_records(l, columns=header)

    # pdb.set_trace()

    round1 = [i for i in set(df.name) if not '-round2' in i]
    round2 = [i for i in set(df.name) if '-round2' in i]

    sm1 = smallest_penalty({k: v for k, v in cand.items() if k in round1})
    sm2 = smallest_penalty({k: v for k, v in cand.items() if k in round2})
    
    selection = df[(df['name'] == sm1) | (df['name'] == sm2)]
    selection.to_csv(outfile, index=False, sep='\t', header=None, quoting=csv.QUOTE_ALL)
    # selection.to_excel(outfile)
    


if __name__ == '__main__':
    validate_primers()

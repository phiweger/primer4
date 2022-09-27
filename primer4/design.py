from collections import defaultdict
from itertools import chain, product
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd
import primer3
from tqdm import tqdm

from primer4.models import PrimerPair
from primer4.utils import log


def design_primers(masked, constraints, params, previous=[]):
    '''
    # 'SEQUENCE_EXCLUDED_REGION'
    # https://primer3.org/manual.html#SEQUENCE_EXCLUDED_REGION
    
    # 'SEQUENCE_INCLUDED_REGION'
    
    # 'PRIMER_PRODUCT_SIZE_RANGE'
    # https://github.com/libnano/primer3-py/issues/18

    No primers? primer3 excludes Ns by default. Relax this and other conditions:

    - https://primer3.org/manual.html#PRIMER_EXPLAIN_FLAG
    - https://primer3.ut.ee/primer3web_help.htm#PRIMER_MAX_NS_ACCEPTED

    Note that the recursive way we implement this here is greedy, i.e. we take
    the best, then mask, then take the next best, ie no "global optimization".
    '''
    # print('Another round!')
    size_range = constraints['size_range']

    # https://primer3.ut.ee/primer3web_help.htm
    # SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=100,50,300,50 ; 900,60,, ; ,,930,100
    only_here = list(chain(*constraints['only_here']))
    # print(only_here)
    # pdb.set_trace()

    if previous:
        x = previous[-1]
        start, end = x.fwd.start, x.fwd.end
        mask_fwd = masked[:start] + ('N' * (end - start)) + masked[end:]

        start, end = x.rev.start, x.rev.end
        mask_rev = masked[:start] + ('N' * (end - start)) + masked[end:]

        masked = ''.join(
            [i if i==j else 'N' for i, j in zip(mask_fwd, mask_rev)])

    spec =  [
        {
            'SEQUENCE_TEMPLATE': masked,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': only_here,
        },
        {
            'PRIMER_NUM_RETURN': params['n_return'],
            'PRIMER_MIN_SIZE': params['size_min'],
            'PRIMER_OPT_SIZE': params['size_opt'],
            'PRIMER_MAX_SIZE': params['size_max'],

            'PRIMER_MAX_NS_ACCEPTED': params['Ns_max'],
    
            'PRIMER_MIN_TM': params['tm_min'],
            'PRIMER_OPT_TM': params['tm_opt'],
            'PRIMER_MAX_TM': params['tm_max'],
            'PRIMER_MIN_GC': params['GC_min'],
            'PRIMER_MAX_GC': params['GC_max'],
    
            'PRIMER_MAX_POLY_X': params['homopolymer_max_len'],
            'PRIMER_MAX_END_GC': params['3prime_max_GC'],
            'PRIMER_MAX_END_STABILITY': params['3prime_stability'],
            'PRIMER_MAX_NS_ACCEPTED': params['Ns_max'],
            
            'PRIMER_PRODUCT_SIZE_RANGE': [size_range],
    
            # defaults, here to be explicit
            'PRIMER_SALT_MONOVALENT': params['salt_monovalent'],
            'PRIMER_SALT_DIVALENT': params['salt_divalent'],
            'PRIMER_DNTP_CONC': params['conc_dNTP'],
            'PRIMER_DNA_CONC': params['conc_DNA'],
        }
    ]

    # https://libnano.github.io/primer3-py/quickstart.html#workflow
    design = primer3.bindings.designPrimers(*spec)
    
    try:
        best = PrimerPair(parse_design(design, 1)[0])
        previous.append(best)
        yield from design_primers(masked, constraints, params, previous)

    except KeyError:
        # No primers found
        yield previous


def parse_design(design, n):
    
    primers = {}
    for i in range(n):
        
        fwd_start, fwd_len = design[f'PRIMER_LEFT_{i}']
        rev_start, rev_len = design[f'PRIMER_RIGHT_{i}']
        
        # c .. candidate
        primers[i] = {
            'fwd': {
                'start': fwd_start,
                'end': fwd_start + fwd_len,
                # 'sanity': template[fwd_start:fwd_start + fwd_len],
                'sequence': design[f'PRIMER_LEFT_{i}_SEQUENCE'],
                'Tm': round(design[f'PRIMER_LEFT_{i}_TM'], 2),
            },
            'rev': {
                # That Python 0-based, end-exclusive indexing thing ...
                'start': rev_start - rev_len + 1,
                'end': rev_start + 1,
                # 'sanity': rc(template[rev_start - rev_len + 1:rev_start + 1]),
                'sequence': design[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                'Tm': round(design[f'PRIMER_RIGHT_{i}_TM'], 2),
            },
            'insert': design[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
            'penalty': round(design[f'PRIMER_PAIR_{i}_PENALTY'], 4),
        }

    return primers


def sort_penalty(primers):
    '''
    It seems primer3 already sorts primers from best to worst, but just to
    make sure.
    '''
    loss = {}
    for k, v in primers.items():
        loss[k] = v['penalty']
    return [k for k, v in sorted(loss.items(), key=lambda x: x[1])]


def check_for_multiple_amplicons(primers, fp_genome, word_size=13, mx_evalue=100, mx_amplicon_len=4000, mx_amplicon_n=1, mn_matches=15, n_cpus=8, mx_blast_hits=10000):
    '''
    Params mostly from ISPCR from UCSC genome browser:

    word_size = 13
    mx_evalue = 10
    mx_amplicon_len = 4000
    mx_amplicon_n = 1  # pseudogene and unwanted amplification check
    mn_matches = 15    # 3' matches

    Testing:

    check_for_multiple_amplicons(
        ('GAGGAAGGCTTTTCGGCATC', 'CTGCGGGAAGCACAGGACAC'),
        'path/to/primer4/data/GRCh37_latest_genomic.fna',
        test=True)
    '''
    tmpdir = tempfile.TemporaryDirectory()
    p = tmpdir.name

    for primer in primers:
        primer.save(f'{p}/{primer.name}.fna')

    fields = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident btop sstrand'
    steps = [
        f'cat {p}/*.fna > {p}/pseudo.fna',
        f'blastn -dust no -word_size {word_size} -evalue {mx_evalue} -outfmt "6 {fields}" -query {p}/pseudo.fna -db {fp_genome} -num_threads {n_cpus} -out {p}/result'
    ]
    '''
    > The “Blast trace-back operations” (BTOP) string describes the alignment produced by BLAST. -- https://www.ncbi.nlm.nih.gov/books/NBK569862

    Examples: 7AG39, 7A-39, 6-G-A41
    '''
    command = '; '.join(steps)
    print(log('Searching for alternative binding sites'))
    log_blast = subprocess.run(command, capture_output=True, shell=True)
    # print(log_blast)

    result = Path(p) / 'result'
    df = pd.read_csv(result, sep='\t', names=fields.split(' '))
    df = df.astype({'btop': str})

    # print(df)
    drop_these = []
    for ix, i in df.iterrows():
        try:
            # Try to parse the number of 3' matches (get the last number)
            # https://stackoverflow.com/questions/5320525/regular-expression-to-match-last-number-in-a-string
            # print(i['btop'])

            matches_3prime = int(re.match(r'.*?(\d+)(?!.*\d)', i['btop']).group(1))
            if matches_3prime < mn_matches:
                drop_these.append(ix)
        except ValueError:
            # Gap or mismatch in last position, so cannot cast string to int
            # invalid literal for int() with base 10: ...
            continue

    # print(f'Dataframe before: {len(df)}')
    df.drop(df.index[drop_these], inplace=True)
    # print(f'Dataframe after: {len(df)}')
    
    unique = set([i.split('.')[0] for i in df['qseqid']])
    
    print(log('Excluding non-unique sites'))
    cnt = defaultdict(int)
    for u in unique:

        fwd_df = df[df['qseqid'] == f'{u}.fwd']
        rev_df = df[df['qseqid'] == f'{u}.rev']

        # If there are too many Blastn hits for a primer pair, the fwd-rev
        # combinations explode thanks to combinatorics. So we make a hard cut
        # and ignore these "abundant" sequences.
        if (len(fwd_df) > mx_blast_hits) or (len(rev_df) > mx_blast_hits):
            cnt[u] = np.Inf
            continue

        fwd = [dict(v) for _, v in fwd_df.iterrows()]
        rev = [dict(v) for _, v in rev_df.iterrows()]
    
        # pr = list(product(fwd, rev))
        # for i, j in tqdm(pr, total=len(pr)):
        for i, j in product(fwd, rev):
            # Same contig?
            if not i['sseqid'] == j['sseqid']:
                continue
    
            # Same chromosome and different strands?
            fwd_strand, fwd_start = i['sstrand'], i['sstart']
            rev_strand, rev_start = j['sstrand'], j['sstart']
            
            # Primers on different strands?
            if not fwd_strand != rev_strand:
                continue
    
            # Orientation right (otherwise polymerase walks in opposite directions)?
            if fwd_strand == 'plus':
                if not fwd_start < rev_start:
                    continue
            else:
                if not fwd_start > rev_start:
                    continue
    
            # Amplicon size?
            if i['sstart'] > j['sstart']:
                amplicon = i['sstart'] - j['send']
            else:
                amplicon = j['sstart'] - i['send']
    
            if amplicon < mx_amplicon_len:
                cnt[u] += 1
                # print(i)
                # print(j)
                # print(amplicon)
                # print(i['qseqid'], i['sseqid'], i['sstart'], i['send'])
                # print(j['qseqid'], j['sseqid'], j['sstart'], j['send'])

    # We should only find one pair for each, which is the amplicon we want.
    results = []
    for primer in primers:
        if cnt[primer.name] <= mx_amplicon_n:
            results.append(primer)
        else:
            print(f'Primer pair {primer.name} does not pass, amplicons calculated: {cnt[primer.name]}')
    # print(results)
    return results


def parse_blast_btop(s):
    '''
    Blast BTOP string .. think sam CIGAR string, but more flexible

    > The “Blast trace-back operations” (BTOP) string describes the alignment
    produced by BLAST. This string is similar to the CIGAR string produced in 
    SAM format, but there are important differences. BTOP is a more flexible 
    format that lists not only the aligned region but also matches and 
    mismatches. BTOP operations consist of 1.) a number with a count of 
    matching letters, 2.) two letters showing a mismatch (e.g., “AG” means A 
    was replaced by G), or 3.) a dash (“-“) and a letter showing a gap. The box 
    below shows a blastn run first with BTOP output and then the same run with 
    the BLAST report showing the alignments.

    -- https://www.ncbi.nlm.nih.gov/books/NBK569862/

    Example: 3GA13, 3GT14, 14, 13-G1TC5, ...

    Usage:

    parse_blast_btop('3GA-13')
    # '...||.............'
    '''
    split = re.findall(r'[A-Za-z]+|\d+|-', s)
    # '3GA-13' > ['3', 'GA', '-', '13']
    result = ''
    for i in split:
        try:
            result += int(i) * '.'
        except ValueError:  # invalid literal for int()
            result += '|'

    return result






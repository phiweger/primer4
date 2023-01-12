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
    # The associated value must be a semicolon-separated list of
    # <left_start>,<left_length>,<right_start>,<right_length>
    # SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=100,50,300,50 ; 900,60,, ; ,,930,100
    only_here = list(chain(*constraints['only_here']))
    # print(only_here)
    # pdb.set_trace()

    # Mask the position of the primer pair from the previous cycle
    if previous:
        x = previous[-1]
        start, end = x.fwd.start, x.fwd.end
        mask_fwd = masked[:start] + ('N' * (end - start)) + masked[end:]

        start, end = x.rev.start, x.rev.end
        mask_rev = masked[:start] + ('N' * (end - start)) + masked[end:]

        masked = ''.join(
            [i if i==j else 'N' for i, j in zip(mask_fwd, mask_rev)])

    #print(only_here)
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
    designs = primer3.bindings.designPrimers(*spec)
    
    try:
        best = PrimerPair(parse_designs(designs, n=1)[0])
        #print(best)
        #from collections import Counter
        #print('Ns:', Counter(masked)['N'])

        if not 'snvs' in constraints:
            raise ValueError('No constraints')
        # DEPRECATED: This will err out if there are no SNVs ie 'snvs': {}
        # if not constraints.get('snvs'):
        #     raise ValueError('No constraints')

        # import pdb; pdb.set_trace()
        d = project_mask_onto_primers(
            best, constraints['snvs'], params['primers']['mn_3prime_matches'])
        
        valid_fwd, dots_fwd, pos_fwd = d['fwd']
        valid_rev, dots_rev, pos_rev = d['rev']
        #print(valid_fwd, dots_fwd, pos_fwd, valid_rev, dots_rev, pos_rev)
        
        if not all([valid_fwd, valid_rev]):
            # We detected an SNV in a primer.
            masked = ''.join(
                ['N' if ix in set(pos_fwd + pos_rev) else i for ix, i in enumerate(masked)])
        else:
            previous.append(best)

        yield from design_primers(masked, constraints, params, previous)

    except KeyError:
        # No primers found
        yield previous


def parse_designs(designs, n):
    
    primers = {}
    for i in range(n):
        
        fwd_start, fwd_len = designs[f'PRIMER_LEFT_{i}']
        rev_start, rev_len = designs[f'PRIMER_RIGHT_{i}']
        
        # c .. candidate
        primers[i] = {
            'fwd': {
                'start': fwd_start,
                'end': fwd_start + fwd_len,
                # 'sanity': template[fwd_start:fwd_start + fwd_len],
                'sequence': designs[f'PRIMER_LEFT_{i}_SEQUENCE'],
                'Tm': round(designs[f'PRIMER_LEFT_{i}_TM'], 2),
            },
            'rev': {
                # That Python 0-based, end-exclusive indexing thing ...
                'start': rev_start - rev_len + 1,
                'end': rev_start + 1,
                # 'sanity': rc(template[rev_start - rev_len + 1:rev_start + 1]),
                'sequence': designs[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                'Tm': round(designs[f'PRIMER_RIGHT_{i}_TM'], 2),
            },
            'insert': designs[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
            'penalty': round(designs[f'PRIMER_PAIR_{i}_PENALTY'], 4),
        }

    return primers


def sort_by_penalty(primers):
    '''
    It seems primer3 already sorts primers from best to worst, but just to
    make sure.
    '''
    loss = {}
    for p in primers:
        loss[p.name] = (p.penalty, p)
    return [v[1] for k, v in sorted(loss.items(), key=lambda x: x[1][0])]


def dereplicate(primers):
    '''
    When we run two design cycles, one allowing SNVs and the other not,
    the name primers can be found independently. Remove those duplicates.
    '''
    duplicate, results = set(), []
    for p1 in primers:
        has_replicate = 0
        d1 = [v for k, v in sorted(p1.data.items())]

        for p2 in primers:
            if p1.name == p2.name:
                continue
            else:
                d2 = [v for k, v in sorted(p2.data.items())]
                if d1 == d2:
                    duplicate.add(tuple(sorted((p1.name, p2.name))))

    # We can only get pairs of duplicates, not 3 duplicates
    if duplicate:
        _, not_those = zip(*duplicate)
    else:
        # No duplicates found
        not_those = set()
    for p in primers:
        if not p.name in not_those:
            # singletons and first in a pair of duplicates
            results.append(p)

    return results


def check_for_multiple_amplicons(primers, fp_genome, params):
    '''
    Params mostly from ISPCR from UCSC genome browser:

    word_size = 13
    mx_evalue = 10
    mx_amplicon_len = 4000
    mx_amplicon_n = 1  # pseudogene and unwanted amplification check
    mn_matches = 15    # 3' matches
    '''
    mx_amplicon_len = params['primers']['mx_amplicon_len']
    mx_amplicon_n = params['primers']['mx_amplicon_n']
    mn_matches = params['primers']['mn_3prime_matches']
    
    word_size = params['blast']['word_size']
    mx_evalue = params['blast']['mx_evalue']
    mx_blast_hits = params['blast']['mx_blast_hits']
    n_cpus = params['blast']['n_cpus']

    blast_index = params['blast']['index']

    tmpdir = tempfile.TemporaryDirectory()
    p = tmpdir.name

    for primer in primers:
        primer.save(f'{p}/{primer.name}.fna')

    fields = 'qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore nident btop sstrand'
    steps = [
        f'cat {p}/*.fna > {p}/pseudo.fna',
        f'blastn -dust no -word_size {word_size} -evalue {mx_evalue} -outfmt "6 {fields}" -query {p}/pseudo.fna -db {blast_index} -num_threads {n_cpus} -out {p}/result'
    ]
    '''
    > The “Blast trace-back operations” (BTOP) string describes the alignment produced by BLAST. -- https://www.ncbi.nlm.nih.gov/books/NBK569862

    Examples: 7AG39, 7A-39, 6-G-A41
    '''
    command = '; '.join(steps)
    print(log(f'Search alternative binding sites for {len(primers)} pair(s)'))
    _ = subprocess.run(command, capture_output=True, shell=True)
    # log_blast = ...
    # print(log_blast)

    result = Path(p) / 'result'
    df = pd.read_csv(result, sep='\t', names=fields.split(' '))
    df = df.astype({'btop': str})
    # print(df.iloc[0])
    # print(f'Before: {len(set([i.split(".")[0] for i in df["qseqid"]]))}')
    # print(len(df))
    # print(df)
    drop_these = []
    aln_profiles = []
    for ix, i in df.iterrows():
        try:
            # Try to parse the number of 3' matches (get the last number)
            # https://stackoverflow.com/questions/5320525/regular-expression-to-match-last-number-in-a-string
            # print(i['btop'])
            profile, split = parse_blast_btop(i['btop'], debug=True)
            aln_profiles.append(profile)
            matches_3prime = split[-1]
            # matches_3prime = int(re.match(r'.*?(\d+)(?!.*\d)', i['btop']).group(1))

            # Is the 3' end of the primer aligned?
            # qend .. End of alignment in query
            # https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
            cond1 = i['qend'] != i['qlen']

            # Are the last <mn_matches> bases matches?
            cond2 = matches_3prime < mn_matches

            # (len(profile) != i['qlen']) .. no softclipping of ends allowed
            #if (matches_3prime < mn_matches) or (len(profile) != i['qlen']):
            
            # 4 debugging
            log_profiles = [
                    matches_3prime,
                    i['qlen'],
                    i['qend'],
                    i['length'],
                    len(split),
                    i['sstrand'],
                    i['btop'],
                    cond1,
                    cond2,
                    profile]

            if cond1 or cond2:
                # TODO: we discard too many? orientation blast/ primer ok?
                # print(log_profiles)
                
                # Primers which bind suboptimally everywhere are thus removed
                drop_these.append(ix)
                
            # else:
            #     print(log_profiles)

        except ValueError:
            # Gap or mismatch in last position, so cannot cast string to int
            # invalid literal for int() with base 10: ...
            continue

    df['aln'] = aln_profiles
    # print(f'Dataframe before: {len(df)}')
    df.drop(df.index[drop_these], inplace=True)
    # print(df)
    # print(f'Dataframe after: {len(df)}')
    # print(f'Remaining: {len(set([i.split(".")[0] for i in df["qseqid"]]))}')
    # print(len(df))
    '''
    qseqid        sseqid  pident  length  ...  nident  btop  sstrand                   aln
    0    139b465b.fwd  NC_000017.10   100.0      20  ...      20    20     plus  ....................
    1    139b465b.fwd  NC_000017.10   100.0      15  ...      15    15     plus       ...............
    '''
    sub = df[['qseqid', 'sseqid', 'sstart', 'send', 'aln']].drop_duplicates()
    lu = {}
    for i in sub.itertuples():
        name, orient = i.qseqid.split('.')
        start, end = i.sstart, i.send
        if start > end:
            start, end = end, start
        lu[(name, orient, i.sseqid, start, end)] = i.aln


    # lu = {tuple(k.split('.')): v for k, v in zip(sub.qseqid, sub.aln)}  
    # Lookup for alignment string:
    # ('139b465b', 'rev'): '...||................'

    # import pdb; pdb.set_trace()

    unique = set([i.split('.')[0] for i in df['qseqid']])
    
    print(log('Exclude non-unique sites'))
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
    # import pdb; pdb.set_trace()
    results = []
    for primer in primers:
        if cnt[primer.name] <= mx_amplicon_n:
            results.append(primer)
        else:
            print(f'Primer pair {primer.name} does not pass, likely too unspecific (> {mx_blast_hits} Blast hits)')
    # print(results)
    return results, lu
    # lu looks like this:
    #  ... ('9434526f', 'rev', 'NT_113943.1', 41544, 41558): '..........|.|..',


def parse_blast_btop(s, debug=False):
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

    parse_blast_btop('3GCA-13')
    # '...||.............'

    Comments:

    The BTOP string will only cover __the aligned fraction__ of the primer! This
    means we actually need to check that the length of the BTOP == len of primer
    or else we miss "soft clippings".

    Example:

    Score = 27.9 bits (14), Expect = 1.6
    Identities = 14/14 (100%)
    Strand = Plus / Minus
                          
    Query: 8        gactgactttctgc 21
                    ||||||||||||||
    Sbjct: 17227357 gactgactttctgc 17227344

    BTOP is in 5' - 3' direction but misses the "ends", which we assume
    non-binding (soft-clipping). See notes.
    '''
    split = re.findall(r'[A-Z][A-Z]|\d+|-[A-Z]|[A-Z]-', s)
    # '3GA-13' > ['3', 'GA', '-', '13']

    cast_split = []
    result = ''

    for i in split:
        try:
            result += int(i) * '.'
            cast_split.append(int(i))
        except ValueError:  # invalid literal for int()
            result += '|'
            cast_split.append(i)
    if debug:
        return result, cast_split
    else:
        return result


def project_mask_onto_primers(primers, mask, mn_3prime_matches=15):
    '''
    Project SNVs onto primer, and decide if compatible with the constraint to
    have a certain number of matches on the 3' end.

    Example:

    project_mask_onto_primers(primers_nomask[0], tmp.mask)
    {
        'fwd': (True, '.....................'),
        'rev': (True, '...||................')
    }
    '''
    d = {}
    if not mask:
        mask = set()

    for x in ['fwd', 'rev']:
        start = primers.data[x]['start']
        end = primers.data[x]['end']
        #print(x, start, end)
        
        z = [('|', i) if i in mask else ('.', i) for i in range(start, end)]
        # import pdb; pdb.set_trace()
        # Dot notation: "." means no SNV, "|" means SNV, e.g. 2nd and last
        # positions of the primer have an SNV: ".|......|"
        dots =  ''.join([i for i, _ in z])
        pos = [j for i, j in z if i == '|']
        assert len(dots) == len(primers.data[x]['sequence']), \
            f'{x}, {dots}, {primers.data[x]["sequence"]}'

        valid = not '|' in dots[-mn_3prime_matches:]
        d[x] = (valid, dots, pos)
    return d





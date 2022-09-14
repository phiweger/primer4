from collections import Counter
import itertools
from io import BytesIO
from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st


def prepare_mock_data_for_vis():
    # Data for plotting
    # https://docs.streamlit.io/library/api-reference/charts/st.pyplot
    return np.random.normal(1, 1, size=300)


class Beauty():
    def __init__(self, data):
        self.data = data
        
    def plot(self, width=10, height=4, dpi=300):
        fig, ax = plt.subplots(figsize=(width, height))
        # fig.set_size_inches(1, 1)
        fig.set_dpi(dpi)
        ax.hist(self.data, bins=30)

        # Oh FFS, there's always something:
        # discuss.streamlit.io/t/cannot-change-matplotlib-figure-size/10295/8
        # TL;DR: st.pyplot(fig) would be ideal but does not scale img correctly
        
        # Workaround for streamlit to respect figure size in matplotlib objects
        buf = BytesIO()
        fig.savefig(buf, format="png")
        # Workaround for streamlit to center the image
        col1, col2, col3 = st.columns([1, 1000, 1])
        with col1:
            st.write('')
        with col2:
            st.image(buf)
        with col3:
            st.write('')    

        return None


def windows(iterable, length=2, overlap=0, truncate=False):
    '''
    Returns a generator of windows of <length> and with an <overlap>.

    Shamelessly stolen from: Python cookbook 2nd edition, chapter 19
    '''

    it = iter(iterable)
    results = list(itertools.islice(it, length))
    
    while len(results) == length:
        yield results
        results = results[length-overlap:]
        results.extend(itertools.islice(it, length-overlap))
    
    if truncate:
        if results and len(results) == length:
            yield results
    else:
        if results:
            yield results


def calculate_gc_content(seq, length=50, overlap=49, truncate=False):
    result = []
    for w in windows(seq.upper(), length, overlap, truncate):
        cnt = Counter(w)
        gc = round((cnt['C'] + cnt['G']) / len(w), 4)
        result.append(gc)

    # Pad last positions w/ 0
    pad = [0] * (len(seq) - len(result))
    _ = result.extend(pad)
    return result


def prepare_data_for_vis(v, tmp, primers, outdir):
    fp = Path(outdir)

    # --- Primers ---
    # bed format
    cnt = 0

    for pair in primers:
        cnt += 1
        result = ''

        for x in ['fwd', 'rev']:
            chrom = tmp.feat.chrom
            start = tmp.invert_relative_pos(pair.data[x]['start'])
            end = tmp.invert_relative_pos(pair.data[x]['end'])
            if not start < end:
                start, end = end, start
            name = uuid4().__str__() + f'::{x}'
            score = cnt
            result += f'{chrom}\t{start}\t{end}\t{name}\t{score}\t.\n'

        with open(fp / f'p{cnt}.bed', 'w+') as out:
            out.write(result)

    # --- Exons ---
    # bed format
    # By default, tmp.region contains all exons.
    features = {}
    for number, i in tmp.exons.items():
        # Exons are already sorted by genomic position from smallest to largest
        start, end = sorted([i.start, i.end])
        features[start] = \
            f'{i.chrom}\t{start}\t{end}\t{number}\t.\t{i.strand}\n'
        bed = ''.join([v for k, v in sorted(features.items())])

    with open(fp / 'exons.bed', 'w+') as out:
        out.write(bed)
    '''
    NC_000017.10    7571739 7573008 11  .   -
    NC_000017.10    7573927 7574033 10  .   -
    NC_000017.10    7576853 7576926 9   .   -
    NC_000017.10    7577019 7577155 8   .   -
    NC_000017.10    7577499 7577608 7   .   -
    NC_000017.10    7578177 7578289 6   .   -
    NC_000017.10    7578371 7578554 5   .   -
    NC_000017.10    7579312 7579590 4   .   -
    NC_000017.10    7579700 7579721 3   .   -
    NC_000017.10    7579839 7579940 2   .   -
    NC_000017.10    7590695 7590808 1   .   -
    '''
    
    # --- SNVs and mask ---
    # bigwig format, 0-based:
    # https://www.biostars.org/p/354635/
    # https://www.biostars.org/p/84686/
    
    # SNVs
    result = ''
    for rel_pos, snvs in tmp.mask_freqs.items():
        gen_pos = tmp.invert_relative_pos(rel_pos)
        for db, freq in snvs:
            # TODO minmax scalar
            freq = np.abs(np.log(freq))
            #freq = 1
            result += f'{tmp.feat.chrom}\t{gen_pos}\t{gen_pos+1}\t{freq}\n'
    
    with open(fp / 'variants.bedgraph', 'w+') as out:
        out.write(result)

    # Mask
    result = ''
    for rel_pos in tmp.mask:
        gen_pos = tmp.invert_relative_pos(rel_pos)
        result += f'{tmp.feat.chrom}\t{gen_pos}\t{gen_pos+1}\t{1}\n'
    
    with open(fp / 'mask.bedgraph', 'w+') as out:
        out.write(result)

    # GC content
    gc = calculate_gc_content(tmp.sequence, 50, 49, truncate=False)
    result = ''
    for rel_pos, v in enumerate(gc):
        gen_pos = tmp.invert_relative_pos(rel_pos)
        result += f'{tmp.feat.chrom}\t{gen_pos}\t{gen_pos+1}\t{v}\n'
    
    with open(fp / 'gc.bedgraph', 'w+') as out:
        out.write(result)

    return None


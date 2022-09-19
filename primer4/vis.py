from collections import Counter
from itertools import islice
from io import BytesIO
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from uuid import uuid4

from jinja2 import Template
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
    results = list(islice(it, length))
    
    while len(results) == length:
        yield results
        results = results[length-overlap:]
        results.extend(islice(it, length-overlap))
    
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


def prepare_data_for_vis(v, tmp, primers):
    '''
    Main fn to prepare, well, data for vis ...
    '''
    max_n_primers = 3

    # --- Primers and query variant ---
    cnt = 0
    tmp_dir = TemporaryDirectory()
    tmp_fp = Path(tmp_dir.name)

    with open(tmp_fp / 'query.bed', 'w+') as out:
        # TODO: this will break,  bc/ c_to_g, what about non-coding?
        # For non-coding variants, we have to apply the offset once the
        # coding position has been mapped to genomic.
        start = tmp.c_to_g[v.start] + v.start_offset
        end = tmp.c_to_g[v.end] + v.end_offset
        # .bed fmt requires that start be smaller then end position
        if start > end:
            start, end = end, start
        # Make non-0 interval if necessary (for plotting to work)
        if start == end:
            end += 1

        line = f'{tmp.feat.chrom}\t{start}\t{end}\t.\t.\t.\n'
        out.write(line)
        # print(line)

    result = ''
    positions = []
    for pair in islice(primers, max_n_primers):
        cnt += 1
        
        for x in ['fwd', 'rev']:
            chrom = tmp.feat.chrom
            
            start = tmp.invert_relative_pos(pair.data[x]['start'])
            end = tmp.invert_relative_pos(pair.data[x]['end'])
            positions.append(start)
            positions.append(end)
            
            if not start < end:
                start, end = end, start
            name = uuid4().__str__() + f'::{x}'
            # score = colors[cnt]
            score = cnt
            result += f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{"+" if x == "fwd" else "-" }\n'
            # print(result)

    with open(tmp_fp / 'primers.bed', 'w+') as out:
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

    with open(tmp_fp / 'exons.bed', 'w+') as out:
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
    
    with open(tmp_fp / 'variants.bedgraph', 'w+') as out:
        out.write(result)

    # Mask
    result = ''
    for rel_pos in tmp.mask:
        gen_pos = tmp.invert_relative_pos(rel_pos)
        result += f'{tmp.feat.chrom}\t{gen_pos}\t{gen_pos+1}\t{1}\n'
    
    with open(tmp_fp / 'mask.bedgraph', 'w+') as out:
        out.write(result)

    # GC content
    gc = calculate_gc_content(tmp.sequence, 50, 49, truncate=False)
    result = ''
    for rel_pos, v in enumerate(gc):
        gen_pos = tmp.invert_relative_pos(rel_pos)
        result += f'{tmp.feat.chrom}\t{gen_pos}\t{gen_pos+1}\t{v}\n'
    
    with open(tmp_fp / 'gc.bedgraph', 'w+') as out:
        out.write(result)


    content_tracks_spec = {
        'query_fp': str(tmp_fp / 'query.bed'),
        'variants_fp': str(tmp_fp / 'variants.bedgraph'),
        'gc_fp': str(tmp_fp / 'gc.bedgraph'),
        'mask_fp': str(tmp_fp / 'mask.bedgraph'),
        'primers_fp': str(tmp_fp / 'primers.bed'),
        'exons_fp': str(tmp_fp / 'exons.bed'),
        'exons_height': 2
    }

    
    with open('tracks.empty.ini', 'r') as file:
        s = file.read()

    empty = Template(s)
    filled = empty.render(**content_tracks_spec)

    with open(tmp_fp / 'tracks.filled.ini', 'w+') as out:
        out.write(filled)

    img_fp = str(tmp_fp / 'img.png')

    subprocess.run([
        'pyGenomeTracks',
        '--tracks', str(tmp_fp / 'tracks.filled.ini'),
        '--region', f'{tmp.feat.chrom}:{np.min(positions)-100}-{np.max(positions)+100}',
        '--outFileName', img_fp
        ])


    from PIL import Image
    image = Image.open(img_fp)

    tmp_dir.cleanup()
    return image

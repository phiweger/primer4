from collections import Counter
from itertools import islice
from io import BytesIO
from pathlib import Path
from PIL import Image
import subprocess
from tempfile import TemporaryDirectory
# from uuid import uuid4

from jinja2 import Template
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from primer4.models import Variant
from primer4.utils import convert_chrom
from primer4.design import project_mask_onto_primers


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
    max_primer_pairs_display = 3
    is_variant = type(v) == Variant

    # --- Primers and query variant ---

    tmp_dir = TemporaryDirectory()
    tmp_fp = Path(tmp_dir.name)

    if type(v) == Variant:
        with open(tmp_fp / 'query.bed', 'w+') as out:
            # TODO: this will break,  bc/ c_to_g, what about non-coding?
            # For non-coding variants, we have to apply the offset once the
            # coding position has been mapped to genomic.
            _, start, end = v.get_genomic_coords(tmp)
            # Make non-0 interval if necessary (for plotting to work)
            if start == end:
                end += 1

            line = f'{tmp.feat.chrom}\t{start}\t{end}\t.\t.\t.\n'
            out.write(line)
            # print(line)

    result = ''
    positions = []
    cnt = 0
    for pair in islice(primers, max_primer_pairs_display):
        # name = uuid4().__str__() + f'::{x}'
        name = cnt
        for x in ['fwd', 'rev']:

            chrom = tmp.feat.chrom
            _, start, end = pair.get_genomic_coords(tmp, x)
            '''
            start = tmp.invert_relative_pos(pair.data[x]['start'])
            end = tmp.invert_relative_pos(pair.data[x]['end'])
            '''
            positions.append(start)
            positions.append(end)

            #if not start < end:
            #    start, end = end, start
            # score = colors[cnt]
            # score = cnt
            score = pair.penalty
            result += f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{"+" if x == "fwd" else "-" }\n'
            # print(result)
        cnt += 1

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

    result = ''
    for rel_pos, snvs in tmp.mask_freqs_filtered.items():
        gen_pos = tmp.invert_relative_pos(rel_pos)
        for db, freq in snvs:
            # TODO minmax scalar
            freq = np.abs(np.log(freq))
            #freq = 1
            result += f'{tmp.feat.chrom}\t{gen_pos}\t{gen_pos+1}\t{freq}\n'
    
    with open(tmp_fp / 'filtered.bedgraph', 'w+') as out:
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
        'filtered_fp': str(tmp_fp / 'filtered.bedgraph'),
        'gc_fp': str(tmp_fp / 'gc.bedgraph'),
        'mask_fp': str(tmp_fp / 'mask.bedgraph'),
        'primers_fp': str(tmp_fp / 'primers.bed'),
        'exons_fp': str(tmp_fp / 'exons.bed'),
    }

    
    with open('tracks.empty.ini', 'r') as file:
        s = file.read()

    empty = Template(s)

    if is_variant:
        comment = ''
    else:
        comment = '#'

    filled = empty.render(**content_tracks_spec, comment=comment)

    with open(tmp_fp / 'tracks.filled.ini', 'w+') as out:
        out.write(filled)

    img_fp = str(tmp_fp / 'img.png')

    subprocess.run([
        'pyGenomeTracks',
        '--tracks', str(tmp_fp / 'tracks.filled.ini'),
        '--region', f'{tmp.feat.chrom}:{np.min(positions)-100}-{np.max(positions)+100}',
        '--outFileName', img_fp,
        '--dpi', '300',
        ],
        stderr=subprocess.DEVNULL)

    try:
        image = Image.open(img_fp)
        tmp_dir.cleanup()
        return image

    except FileNotFoundError:
        tmp_dir.cleanup()
        st.error('No image generated, maybe "tracks.ini" contains errors?')
        st.stop()


def primers_to_df(primers, tmp, qry, aln):
    pd.set_option("display.precision", 2)

    l = []
    for pair in primers:
        _, fwd_start, fwd_end = pair.get_genomic_coords(tmp, 'fwd')
        _, rev_start, rev_end = pair.get_genomic_coords(tmp, 'rev')

        # TODO: I already did this in the recursion that designs the primers;
        # however, for now I accept the slight code duplication. Could add
        # this to the PrimerPair object or keep here so it is explicit.
        d = project_mask_onto_primers(pair, tmp.mask)
        _, dots_fwd, pos_fwd = d['fwd']
        _, dots_rev, pos_rev = d['rev']

        row = [
            pair.name,
            pair.penalty,
            pair.get_amplicon_len(),
            len(pair.data['fwd']['sequence']),
            len(pair.data['rev']['sequence']),
            pair.get_gc('fwd'),
            pair.get_gc('rev'),
            pair.data['fwd']['Tm'],
            pair.data['rev']['Tm'],
            pair.data['fwd']['sequence'],
            pair.data['rev']['sequence'],
            tmp.tx,
            tmp.feat.attributes.get('gene')[0],
            convert_chrom(tmp.feat.chrom),
            fwd_start,
            fwd_end,
            rev_start,
            rev_end,
            *[f'c.{i}' for i in pair.get_coding_coords(tmp, 'fwd')[1:]],
            *[f'c.{i}' for i in pair.get_coding_coords(tmp, 'rev')[1:]],
            qry,
            dots_fwd,
            dots_rev,
            #aln[(pair.name, 'fwd', tmp.feat.chrom, fwd_start, fwd_end)],
            #aln[(pair.name, 'rev', tmp.feat.chrom, rev_start, rev_end)],
        ]
        l.append(row)

    # Any primers found?
    if not l:
        return pd.DataFrame()

    else:
        # TODO: Add mismatches and poistion of mm from 3' end
        df = pd.DataFrame(l)
        df.columns = 'name,penalty,amplicon,fwd len,rev len,fwd GC,rev GC,fwd Tm,rev Tm,fwd 5>3,rev 5>3,transcript,gene,chrom,fwd start,fwd end,rev start,rev end,fwd c. start, fwd c. end, rev c. start, rev c. end,query,aln fwd 5>3,aln rev 5>3'.split(',')    
        # Sort df; in case of qPCR we look first left then right of exon so we
        # get two independent sets of primers, ie df is not ordered in this case
        df = df.sort_values('penalty')
        return df


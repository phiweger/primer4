from io import BytesIO

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


def prepare_data_for_vis(v, tmp, primers, prefix):
    
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

    with open(f'{prefix}.bed', 'w+') as out:
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
    
    # --- SNVs ---
    # bigwig format, 0-based:
    # https://www.biostars.org/p/354635/
    # https://www.biostars.org/p/84686/
    result = ''
    for rel_pos, snvs in tmp.mask_freqs.items():
        gen_pos = tmp.invert_relative_pos(rel_pos)
        for db, freq in snvs:
            freq = np.log(freq)
            result += f'{tmp.feat.chrom}\t{gen_pos}\t{gen_pos+1}\t{freq}\n'
    
    with open(f'{prefix}.bw', 'w+') as out:
        out.write(result)

    return None


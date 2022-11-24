'''
streamlit run primer4/stream.py -- --config config.json
pytest --config config.json
'''


import argparse
# from collections import Counter
import json
from pathlib import Path
import pdb
import os

import click
from cdot.hgvs.dataproviders import JSONDataProvider
import gffutils
import hgvs
import pandas as pd
from pyfaidx import Fasta
import streamlit as st

from primer4.models import Variant, ExonDelta, SingleExon, ExonSpread, Template
from primer4.design import (
    design_primers,
    check_for_multiple_amplicons,
    sort_by_penalty,
    dereplicate,
    )
from primer4.hacks import download_button
from primer4.utils import (
    log,
    mask_sequence,
    reconstruct_mrna,
    sync_tx_with_feature_db,
    )
from primer4.vis import prepare_data_for_vis, primers_to_df
# from primer4.vis import Beauty, prepare_mock_data_for_vis
from primer4.warnings import warn


def gimme_some_primers(method, code, fp_genome, genome, hdp, db, vardbs, params, max_variation, blind_search):
    '''
    gimme_some_primers('sanger', 'NM_000546.6:c.215C>G', ...)
    gimme_some_primers('qpcr', ('NM_001145408.2', 6), ...)
    gimme_some_primers('mrna', ('NM_000546.6', 6, 7), ...)
    '''
    if method == 'sanger':

        # This would otherwise fail later on HGVS parsing error
        if len(code) > 1:
            raise ValueError('Wrong query syntax, should be something like "NM_000546.6:c.215C>G"')
        
        try:
            v = Variant(code[0], hdp, db)
        # TODO: Catch this error in the main fn so we don't have st.x() fn
        # all over the code? We could let all errors propagate up to there ...
        except hgvs.exceptions.HGVSParseError:
            st.warning('Variant could not be parsed, check syntax (spaces?).')
            st.stop()

    elif method == 'qpcr':
        if len(code) != 2:
            raise ValueError(f'Wrong query syntax, should be like "NM_000546.6::4"')
        v = SingleExon(code[0], int(code[1]))
    
    elif method == 'mrna':
        if len(code) != 3:
            raise ValueError(f'Wrong query syntax, should be like "NM_000546.6::5::6"')   
        v = ExonSpread(code[0], int(code[1]), int(code[2]))
    
    else:
        raise ValueError('Method is not implemented, exit.')

    tmp = Template(v, db)
    
    # Annotate templete
    tmp.load_variation_freqs_(vardbs, params)
    tmp.load_variation_(max_variation)
    tmp.get_sequence_(genome)

    # Mask and get primers
    if method == 'sanger':
        
        constraints = tmp.apply(method, db, params)
        masked = mask_sequence(tmp.sequence, tmp.mask)
        
        # We can run primer4 in two ways: First mark SNVs, then search primers
        # "between" them OR first search primers, and then filter them if any
        # SNVs are located in sensitive areas such as 3' binding sites.
        print(log('Run: SNVs first'))
        constraints['snvs'] = tmp.mask
        primers = [p for p in next(
            design_primers(masked, constraints, params, []))]
        print(log(f'Found {len(primers)} primers'))
        #print(constraints)
        
        if blind_search:
            nomask = mask_sequence(tmp.sequence, set())
            #constraints['snvs'] = tmp.mask
            #import pdb
            #pdb.set_trace()
            print(log('Run: Design first'))
            primers_nomask = [p for p in next(
                design_primers(nomask, constraints, params, []))]
            # print(constraints)
            print(log(f'Found {len(primers_nomask)} more primers'))
            primers = primers + primers_nomask
            print(log(f'Found {len(primers)} in total'))


    elif method == 'qpcr':
        masked = mask_sequence(tmp.sequence, tmp.mask)
        
        print(log('Run: SNVs first'))
        primers = []
        all_constraints = tmp.apply('qpcr', db, params)
        for constraints in all_constraints:
            # print(constraints)
            constraints['snvs'] = tmp.mask
            x = [p for p in next(design_primers(masked, constraints, params, []))]
            primers.extend(x)
        print(log(f'Found {len(primers)}'))

        if blind_search:
            nomask = mask_sequence(tmp.sequence, set())
            #constraints['snvs'] = tmp.mask
            #import pdb
            #pdb.set_trace()
            print(log('Run: Design first'))
            
            for constraints in all_constraints:
                primers_nomask = [p for p in next(
                    design_primers(nomask, constraints, params, []))]
                primers.extend(primers_nomask)
            print(f'Found {len(primers)} in total')


    elif method == 'mrna':
        # tmp.mrna = reconstruct_mrna(tmp, db, genome)
        
        # Mask all exons but the ones we want to find primers in
        tmp.mrna = reconstruct_mrna(tmp, db)
        # Contains: mrna_mask, exons, offset
        
        # Mask SNVs
        masked = mask_sequence(tmp.sequence, tmp.mask)
        # Now merge
        mrna_mask, _, offset = tmp.mrna
        assert len(mrna_mask) == len(masked)
        masked = ''.join([j if j=='N' else i for i, j in zip(mrna_mask, masked)])

        constraints = tmp.apply('mrna', db, params)
        constraints['snvs'] = tmp.mask

        print(log('Run: SNVs first'))
        primers = [p for p in next(design_primers(masked, constraints, params, []))]
        # import pdb
        # pdb.set_trace()

        if blind_search:
            nomask = mask_sequence(tmp.sequence, set())
            nomask = ''.join([j if j=='N' else i for i, j in zip(mrna_mask, nomask)])

            print(log('Run: Design first'))
            primers_nomask = [p for p in next(design_primers(nomask, constraints, params, []))]
            # print(constraints)
            print(log(f'Found {len(primers_nomask)} more primers'))
            primers = primers + primers_nomask
            print(log(f'Found {len(primers)} in total'))
        
        # Add offset, inplace operation
        for p in primers:
            p.offset = offset


    else:
        raise ValueError('Method is not implemented, exit.')
    
    
    if blind_search:
        print(log('Dereplicating primers from two search cycles (with and without SNVs)'))
        primers = dereplicate(primers)

    mx = params['primers']['check_max_num_candidates']
    primers = sort_by_penalty(primers)[:mx]

    results, aln = check_for_multiple_amplicons(primers, fp_genome, params)
    # Until now, we have only checked the alignment of primers to the
    # reference genome -- any "variants" are really mapping mismatches.
    # In the case of designing primers while ignoring SNVs, we need to
    # add those SNVs back in.

    #import pdb
    #pdb.set_trace()
    
    msg = f'{len(results)}/{len(primers)} primer pairs pass all filters'
    print(log(msg))
    st.write(msg)
    return results, tmp, aln


def load_chromosome_names(fn):
    # https://stackoverflow.com/questions/1270951/how-to-refer-to-relative-paths-of-resources-when-working-with-a-code-repository
    # fn = Path(__file__).parents[1] / 'chrom_names_hg38.csv'
    chrom_names = {}
    with open(fn, 'r') as file:
        for line in file:
            k, v = line.strip().split(',')
            chrom_names[k] = v
    return chrom_names


# https://docs.streamlit.io/knowledge-base/using-streamlit/caching-issues
# https://discuss.streamlit.io/t/unhashabletype-cannot-hash-object-of-type-thread-local/1917
@st.cache(allow_output_mutation=True)
def housekeeping(params):
    print(log('Housekeeping ...'))

    fp_genome = params['data']['reference']
    fp_coords = params['data']['coordinates']
    vardbs = params['data']['variation']

    # We load the annotation data later (see comment in main fn) but check
    # existance here.
    fp_annotation = params['data']['annotation']

    x = [i for i in vardbs.values()]
    for fp in [fp_coords, fp_annotation, fp_genome] + x:
        p = Path(fp)
        assert p.exists(), f'Path {fp} does not exist'

    print(log('Loading reference genome sequence'))
    genome = Fasta(fp_genome)
    print(log('Loading transcript coordinate mappings'))
    hdp = JSONDataProvider([fp_coords])

    return genome, hdp, vardbs


# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--config', required=True, type=str)
args = parser.parse_args()


def main(fp_config):
    st.set_page_config(layout='centered')

    with open(fp_config, 'r') as file:
        params = json.load(file)
    
    chrom_names = load_chromosome_names(params['data']['chrom_names'])
    _ = params.update({'cn': chrom_names})  # inplace, this feels dirty

    # TODO:
    # https://github.com/phiweger/primer4/issues/2
    # Why does it trigger reload in housekeeping fn but not here?
    # Something w/ session state?
    # https://stackoverflow.com/questions/64698788/streamlit-how-to-store-value-of-variable-in-cache
    # https://github.com/biocommons/biocommons.seqrepo
    x = params['data']['sequences']
    if Path(x).exists():
        print(log('Will use local transcript sequence data'))
        # TODO: It's this line:
        os.environ['HGVS_SEQREPO_DIR'] = x
    else:
        # If env var is not set, hgvs library will default to API usage, ie
        # internet connection is needed.
        print(log('Will use API to obtain sequence data'))

    genome, hdp, vardbs = housekeeping(params)

    # Why is load annotation data here, not in housekeeping fn?
    # Cannot be opened by housekeeping bc/ second iteration will cause:
    # sqlite3.ProgrammingError: SQLite objects created in a thread can only be used in that same thread. The object was created in thread id 123145481936896 and this is thread id 123145532841984.
    db = gffutils.FeatureDB(params['data']['annotation'], keep_order=True)


    st.markdown(
        r'''
        ## ❤️ Primer4U ❤️

        Example queries:

        ```bash
        # Sanger; HGVS syntax
        NM_000546.6:c.215C>G
        NM_005585.4:c.1419dup
        # qPCR; eg "::3" means we target exon 3
        NM_000546.6::4
        # mRNA; anchor primers in two exons
        NM_000546.6::5::7
        ```
        '''
        )

    # The menu
    # https://docs.streamlit.io/library/api-reference/widgets
    # Form
    # https://docs.streamlit.io/library/api-reference/control-flow/st.form
    with st.form('my_form'):
        
        order = st.text_input('Query', '')
        # Based on the selected method (PCR, ...) we'd like to update the
        # default values for the fields below. In theory, this is possible:
        # https://discuss.streamlit.io/t/circular-connection-of-slider-and-text-input/11015/4
        # BUT. Not inside a form:
        # streamlit.errors.StreamlitAPIException: With forms, callbacks can 
        # only be defined on the `st.form_submit_button`. Defining callbacks on 
        # other widgets inside a form is not allowed.
        col1, col2, col3, col4 = st.columns([1, 1, 1, 1])
        with col1:
            method = st.selectbox('Method', ('Sanger', 'qPCR', 'mRNA'))
            # on_change=foo
            # StreamlitAPIException: With forms, callbacks can only be defined 
            # on the st.form_submit_button. Defining callbacks on other widgets 
            # inside a form is not allowed.
            method = method.lower()

        with col2:
            amplicon_len_min = st.number_input('min length [bp]', value=350)
        with col3:
            amplicon_len_max = st.number_input('max length [bp]', value=600)
        with col4:
            max_variation = st.number_input('Allele frequency [%]', min_value=0., max_value=100., value=0., step=0.01, format='%.2f') / 100
        
        # Row 2
        # https://discuss.streamlit.io/t/how-to-have-2-rows-of-columns-using-st-beta-columns/11699/2
        with col1:
            blind_search = st.checkbox('Allow SNVs', value=True)

        # Every form must have a submit button.
        submitted = st.form_submit_button(
            'Run',
            on_click=warn(method, params, amplicon_len_min, amplicon_len_max))
            # warn() will replace params inplace (eg amplicon len)

        if submitted:
            if not order:
                st.write('Please provide a query')
                return None
            else:
                code = order.split('::')  # case: NM_000546.6::3

                # Replace transcript version if necessary
                tx = code[0].split(':')[0]  # case : NM_005585.4:c.1419dup
                
                print(log('Sync transcript'))
                used_tx = sync_tx_with_feature_db(tx, db)
                if used_tx != tx:
                    st.warning(f'Used trancript {used_tx}')
                    code = [code[0].replace(tx, used_tx)] + code[1:]

                print(log('Primer design'))
                with st.spinner(text='Design in progress ...'):
                    primers, tmp, aln = gimme_some_primers(
                        method,
                        code,
                        params['data']['reference'],
                        genome,
                        hdp,
                        db,
                        vardbs,
                        params,
                        max_variation,
                        blind_search)

        else:
            return None
    # What's in "primers"?
    #import pdb
    #pdb.set_trace()
    # if primers: dir(primers[0])
    # dir(tmp)
    # [... 'data', 'fwd', 'insert', 'name', 'penalty', 'rev', 'save', 'to_c', 
    # 'to_g']
    # primers[0].to_g(tmp)
    # [7579197, 7579217, 7579679, 7579699]
    # TODO: coding coords

    # TODO def prepare_df() in vis
    
    df = primers_to_df(primers, tmp, order, params)
    if df.empty:
        st.write('No primers found under the provided constrains. Relax (the constraints)!')

    else:
        # Plot something
        # data = prepare_mock_data_for_vis()
        #_ = Beauty(data).plot()
        
        # "image" is generated like so:
        # from PIL import Image
        # image = Image.open(fp)
        with st.spinner(text='Preparing results ...'):
            image = prepare_data_for_vis(tmp.data, tmp, primers)
    
            # Center image
            col1, col2, col3 = st.columns([1, 1000, 1])
            st.text('\n')
            with col1:
                st.write('')
            with col2:
                st.image(image)
            with col3:
                st.write('')    
            #from PIL import Image
            #image = Image.open(img_fp)
            # st.image(image)
            #img_fp.cleanup()
    
    
            # https://docs.streamlit.io/library/api-reference/data/st.dataframe
            # st.table(df)
            st.text('\n')
            st.dataframe(
                df.style.format(
                    {
                        'penalty': '{:.2f}',
                        'fwd GC': '{:.2f}',
                        'rev GC': '{:.2f}',
                        'fwd Tm': '{:.2f}',
                        'rev Tm': '{:.2f}',
                    }))
        
        # TODO: Style columns?
        # https://stackoverflow.com/questions/41654949/pandas-style-function-to-highlight-specific-columns
        # st.dataframe(df.style.highlight_max(axis=1))
        
        # Download
        # https://docs.streamlit.io/knowledge-base/using-streamlit/how-download-pandas-dataframe-csv
        # https://docs.streamlit.io/knowledge-base/using-streamlit/how-download-file-streamlit
        @st.cache
        def convert_df(df):
            return df.to_csv(index=False, quoting=None).encode('utf-8')

        #import pdb
        #pdb.set_trace()
        csv = convert_df(df)

        # This reloads the entire page, see "hacks.py"
        # st.download_button(
        #     "Download",
        #     csv,
        #     "file.csv",
        #     "text/csv",
        #     key='download-csv'
        #     )

        #import pandas as pd
        
        download_button_str = download_button(csv, 'primers.csv', f'Download')
        st.markdown(download_button_str, unsafe_allow_html=True)

    return None


if __name__ == '__main__':
    main(args.config)


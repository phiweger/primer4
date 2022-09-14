'''
streamlit run primer4/stream.py -- --config config.json
pytest --config config.json
'''


import argparse
import json
from pathlib import Path
import pdb
import os

import click
from cdot.hgvs.dataproviders import JSONDataProvider
import gffutils
import pandas as pd
from pyfaidx import Fasta
import streamlit as st

from primer4.models import Variant, ExonDelta, SingleExon, ExonSpread, Template
from primer4.design import design_primers, check_for_multiple_amplicons
from primer4.utils import mask_sequence, reconstruct_mrna, log
from primer4.vis import Beauty, prepare_mock_data_for_vis
from primer4.warnings import warn


def gimme_some_primers(method, code, fp_genome, genome, hdp, db, vardbs, params, max_variation):
    '''
    gimme_some_primers('sanger', 'NM_000546.6:c.215C>G', ...)
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
    
    # Annotate templete
    tmp.load_variation_freqs_(vardbs)
    tmp.load_variation_(max_variation)
    tmp.get_sequence_(genome)

    # Mask and get primers
    if method == 'sanger':
        masked = mask_sequence(tmp.sequence, tmp.mask)
        constraints = tmp.apply(method, db, params)
        primers = [p for p in next(design_primers(masked, constraints, params, []))]

    elif method == 'qpcr':
        masked = mask_sequence(tmp.sequence, tmp.mask)
        
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
    return results, tmp


# https://docs.streamlit.io/knowledge-base/using-streamlit/caching-issues
# https://discuss.streamlit.io/t/unhashabletype-cannot-hash-object-of-type-thread-local/1917
@st.cache(allow_output_mutation=True)
def housekeeping(fp_config):
    print('Housekeeping ...')

    with open(fp_config, 'r') as file:
        params = json.load(file)

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

    return genome, hdp, params, vardbs


# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--config', required=True, type=str)
args = parser.parse_args()


def main(fp_config):
    genome, hdp, params, vardbs = housekeeping(fp_config)

    # Why is load annotation data here, not in housekeeping fn?
    # Cannot be opened by housekeeping bc/ second iteration will cause:
    # sqlite3.ProgrammingError: SQLite objects created in a thread can only be used in that same thread. The object was created in thread id 123145481936896 and this is thread id 123145532841984.
    db = gffutils.FeatureDB(params['data']['annotation'], keep_order=True)

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


    st.markdown(
        r'''
        ## Primer4

        Example queries:

        ```bash
        # Sanger; HGVS syntax
        NM_000546.6:c.215C>G
        # mRNA; eg "::3" means we target exon 3
        NM_000546.6::3
        # qPCR; anchor primers in two exons
        NM_000546.6::5::6 
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
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            method = st.selectbox('Method', ('Sanger', 'qPCR', 'mRNA'))
            method = method.lower()
        with col2:
            amplicon_len_min = st.number_input('min length [bp]', value=250)
        with col3:
            amplicon_len_max = st.number_input('max length [bp]', value=600)
        with col4:
            max_variation = st.number_input('Allele frequency [%]', min_value=0., max_value=100., value=0., step=0.01) / 100
    
        # Every form must have a submit button.
        submitted = st.form_submit_button(
            'Run',
            on_click=warn(method, params, amplicon_len_min, amplicon_len_max))
        
        if submitted:
            if not order:
                st.write('Please provide a query')
                return None
            else:
                code = order.split('::')
                # This would otherwise fail later on HGVS parsing error
                if len(code) > 1 and method == 'sanger':
                    raise ValueError('Wrong query syntax, should be something like "NM_000546.6:c.215C>G"')
    
                primers, tmp = gimme_some_primers(
                    method,
                    code,
                    params['data']['reference'],
                    genome,
                    hdp,
                    db,
                    vardbs,
                    params,
                    max_variation)
                
                st.write('Done.')
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

    l = []
    for pair in primers:
        # with open(pout / f'{pair.name}.json', 'w+') as out: 
        #     json.dump(pair.data, out, indent=4, sort_keys=True)
        # print(pair.name)
        # print(pair.data)

        l.append(f"{order},{pair.name},{pair.penalty},{pair.data['fwd']['sequence']},{pair.data['fwd']['Tm']},{pair.data['rev']['sequence']},{pair.data['rev']['Tm']}\n")

    # Any primers found?
    if not l:
        st.write('No primers found under the provided constrains. Relax! (the constraints)')
    
    else:
        # Display dataframe
        df = pd.DataFrame([i.split(',') for i in l])
        df.columns = 'order name penalty fwd fwd_tm rev rev_tm'.split(' ')    
        # https://docs.streamlit.io/library/api-reference/data/st.dataframe
        # st.table(df)
        st.dataframe(df)

        # Plot something
        data = prepare_mock_data_for_vis()
        _ = Beauty(data).plot()

        # Download
        # https://docs.streamlit.io/knowledge-base/using-streamlit/how-download-pandas-dataframe-csv
        # https://docs.streamlit.io/knowledge-base/using-streamlit/how-download-file-streamlit
        @st.cache
        def convert_df(df):
            return df.to_csv().encode('utf-8')

        csv = convert_df(df)
        st.download_button(
            "Download",
            csv,
            "file.csv",
            "text/csv",
            key='download-csv'
            )

    return None


if __name__ == '__main__':
    main(args.config)


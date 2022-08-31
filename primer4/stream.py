import json
from pathlib import Path
import pdb

import click
from cdot.hgvs.dataproviders import JSONDataProvider
import gffutils
from pyfaidx import Fasta
import streamlit as st

from primer4.models import Variant, ExonDelta, SingleExon, ExonSpread, Template
from primer4.design import design_primers, check_for_multiple_amplicons
from primer4.utils import mask_sequence, reconstruct_mrna


def gimme_some_primers(method, code, fp_genome, genome, hdp, db, vardbs, params):
    '''
    fp_data = '/Users/phi/Dropbox/repos/primer4/data'
    fp_config = '/Users/phi/Dropbox/repos/primer4/config.json'

    gimme_some_primers('sanger', 'NM_000546.6:c.215C>G', fp_data, fp_config)
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
    tmp.load_variation_(vardbs)
    
    if method == 'sanger':
        masked = mask_sequence(tmp.get_sequence(genome), tmp.mask)
        constraints = tmp.apply(method, db, params)
        primers = [p for p in next(design_primers(masked, constraints, params, []))]

    elif method == 'qpcr':
        masked = mask_sequence(tmp.get_sequence(genome), tmp.mask)
        
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
    return results



import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-i', '--order', required=False)
parser.add_argument('-o', '--outdir', default='results')
parser.add_argument('-d', '--fp-data')
parser.add_argument('-p', '--fp-config')
args = parser.parse_args()


order = args.order
outdir = args.outdir
fp_data = args.fp_data
fp_config = args.fp_config


# https://docs.streamlit.io/knowledge-base/using-streamlit/caching-issues
# https://discuss.streamlit.io/t/unhashabletype-cannot-hash-object-of-type-thread-local/1917
@st.cache(allow_output_mutation=True)
def housekeeping(fp_data, fp_config):
    print('Housekeeping ...')
    p = Path(fp_data)
    assert p.exists(), 'Data path does not exist, exit.'
    
    fp_genome = str(p / 'GRCh37_latest_genomic.fna')
    fp_coords = str(p / 'cdot-0.2.1.refseq.grch37_grch38.json.gz')
    fp_annotation = str(p / 'hg19-p13_annotation_bak.db')
    # fp_annotation = f'{fp_data}/hg19-p13_annotation.db'
    
    fp_snvs_1 = str(p / 'GRCh37_latest_dbSNP_all.vcf.gz')
    fp_snvs_2 = str(p / 'ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz')
    fp_snvs_3 = str(p / 'ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.gz')

    genome = Fasta(fp_genome)
    hdp = JSONDataProvider([fp_coords])
    
    # TODO:
    # os.environ['HGVS_SEQREPO_DIR'] = '.../seqrepo/2021-01-29/'
    
    with open(fp_config, 'r') as file:
        params = json.load(file)
    
    vardbs = {
        'dbSNP': fp_snvs_1,
        '1000Genomes': fp_snvs_2,
        'ESP': fp_snvs_3
        }
    return fp_genome, genome, hdp, params, vardbs


def main(order, outdir, fp_data, fp_config):
    fp_genome, genome, hdp, params, vardbs = housekeeping(fp_data, fp_config)

    # Why is "fp_annotation" here?
    # Cannot be opened by housekeeping bc/ second iteration will cause:
    # sqlite3.ProgrammingError: SQLite objects created in a thread can only be used in that same thread. The object was created in thread id 123145481936896 and this is thread id 123145532841984.
    p = Path(fp_data)
    assert p.exists(), 'Data path does not exist, exit.'
    fp_annotation = str(p / 'hg19-p13_annotation_bak.db')
    db = gffutils.FeatureDB(fp_annotation, keep_order=True)

    pout = Path(outdir)
    pout.mkdir(exist_ok=True)

    # The menu
    # https://docs.streamlit.io/library/api-reference/widgets

    
    # Form
    # https://docs.streamlit.io/library/api-reference/control-flow/st.form

    with st.form("my_form"):
        method = st.selectbox(
            'What method?',
            ('Sanger', 'qPCR', 'mRNA'))
        order = st.text_input('Variant', '')

        st.markdown(
            r'''
            Examples:

            NM_000546.6:c.215C>G
            NM_000546.6::3 (note the double "::")
            ''')

        mask_variants = st.checkbox('Mask variants', value=True)
    
        # Every form must have a submit button.
        submitted = st.form_submit_button('Gimme some primers!')
        
        if submitted:
            if not order:
                return None
            else:
                code = order.split('::')
                method = method.lower()
                st.write(method, code)
    
                primers = gimme_some_primers(method, code, fp_genome, genome, hdp, db, vardbs, params)
                st.write('Done.')
        else:
            return None

    # method, gene, code = order.strip().split('::')

    l = []
    for pair in primers:
        # st.write(pair, pair.data['fwd'], pair.data['rev'])
        # if code[0] == 'NM_006087.4:c.745G>A' and pair.penalty < 1:
            # pdb.set_trace()

        with open(pout / f'{pair.name}.json', 'w+') as out: 
            json.dump(pair.data, out, indent=4, sort_keys=True)
        # print(pair.name)
        # print(pair.data)

        l.append(f"{order},{pair.name},{pair.penalty},{pair.data['fwd']['sequence']},{pair.data['fwd']['Tm']},{pair.data['rev']['sequence']},{pair.data['rev']['Tm']}\n")

    # with open(f'{outdir}.csv', 'w+') as out:
    #     for i in l:
    #         out.write(i)
    

    # TODO: st.table
    # https://docs.streamlit.io/library/api-reference/data

    
    # https://docs.streamlit.io/knowledge-base/using-streamlit/how-download-file-streamlit
    st.download_button('Anneal me', ''.join(l), 'primers.csv')

    return None


if __name__ == '__main__':
    main(order, outdir, fp_data, fp_config)

import json
from pathlib import Path

import pytest

from cdot.hgvs.dataproviders import JSONDataProvider
import gffutils
from pyfaidx import Fasta
from screed import rc

from primer4.models import Variant, Template



# default_data_path = Path('data')

# https://stackoverflow.com/questions/51703577/can-we-pass-html-log-path-in-pytest-from-a-variable-in-script
# @pytest.hookimpl(tryfirst=True)
# def pytest_configure(config):
#     if not config.option.data:
#         config.option.data = default_data_path



# @pytest.fixture
# def setup_option(request):
#     return request.config.getoption("--data")



def test_that(setup_option):
    # print("setup_option: %s" % setup_option)
    print(request)



code = 'NM_000546.6:c.215C>G'

fp_data = '/Users/phi/Dropbox/repos/primer4/data'

p = Path(fp_data)
assert p.exists(), 'Data path does not exist, exit.'

fp_genome = str(p / 'GRCh37_latest_genomic.fna')
fp_coords = str(p / 'cdot-0.2.1.refseq.grch37_grch38.json.gz')
fp_annotation = str(p / 'hg19-p13_annotation_bak.db')
# fp_annotation = f'{fp_data}/hg19-p13_annotation.db'

fp_snvs_1 = str(p / 'GRCh37_latest_dbSNP_all.vcf.gz')
fp_snvs_2 = str(p / 'ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz')
fp_snvs_3 = str(p / 'ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.gz')


hdp = JSONDataProvider([str(fp_coords)])
db = gffutils.FeatureDB(fp_annotation, keep_order=True)
vardbs = {
    'dbSNP': fp_snvs_1,
    '1000Genomes': fp_snvs_2,
    'ESP': fp_snvs_3
    }


v = Variant(code, hdp, db)
tmp = Template(v, db)
tmp.load_variation_(vardbs)

# tmp.get_sequence(Fasta(fp_genome))


# TODO: make SNV allel freq available, other filters?


boundaries = []
for i in db.children(tmp.feat, featuretype='exon', order_by='start'):
    # exon-NM_000546.6-1
    num = int(i.id.split('-')[-1])
    boundaries.append({
        'number': num,
        'start': i.start,
        'end': i.end,
        'strand': i.strand,
        })


vis = {}
vis['transcript'] = tmp.get_sequence(Fasta(fp_genome)).upper()
vis['variation'] = tmp.mask_freqs
vis['chrom'] = tmp.feat.chrom
vis['start'], vis['end'] = tmp.start, tmp.end  # pythonic boundaries!


# Variant
vis['mutation'] = {
    'start': v.start,
    'end': v.end,
    'start_offset': v.start_offset,
    'end_offset': v.end_offset,
}


# Mock primers
primers = []

fwd = vis['transcript'][430:451]
rev = rc(vis['transcript'][500:521])
p1 = {
    'name': 'foo',
    'fwd': {
        'sequence': fwd,
        'start': 430,
        'end': 451-1
        },
    'rev': {
        'sequence': rev,
        'start': 521-1,
        'end': 500,
        }
    }


fwd = vis['transcript'][408:429]
rev = rc(vis['transcript'][550:571])
p2 = {
    'name': 'bar',
    'fwd': {
        'sequence': fwd,
        'start': 408,
        'end': 429-1
        },
    'rev': {
        'sequence': rev,
        'start': 550-1,
        'end': 571,
        }
    }



vis['primers'] = [p1, p2]
vis['exon_boundaries'] = boundaries

with open('4vis.json', 'w+') as out:
    json.dump(vis, out, indent=4)

'''
TODO:

there's an off by one error

    "variation": {
        "-1": [
            [
                "dbSNP",
                1.511e-05
            ]
        ],
        "1": [
            [
                "dbSNP",
                0.0
            ]
        ],

utils.py line 559
'''





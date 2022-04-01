import json
import pytest

from cdot.hgvs.dataproviders import JSONDataProvider
from gffutils import FeatureDB

from primer4.models import Variant, ExonDelta


with open('config.json', 'r') as file:
    params = json.load(file)

hdp = JSONDataProvider([f'data/{params["coordinates"]}'])
db = FeatureDB(f'data/{params["annotation"]}', keep_order=True)


# https://towardsdatascience.com/pytest-for-data-scientists-2990319e55e6
signature = 'code, coords, feature_db, coding_start, chrom, genomic_start'
testdata = [
    ('NM_000546.6:c.215C>G', hdp, db, 215, 'NC_000017.10', 7579472),
    ('NM_000546.6:c.672+3C>G', hdp, db, 672, 'NC_000017.10', 7578174),
    ('NM_015015.3:c.2441+1G>A', hdp, db, 2441, 'NC_000019.9', 5137688)
    ]


@pytest.mark.parametrize(signature, testdata)
def test_variant(code, coords, feature_db, coding_start, chrom, genomic_start):
    v = Variant(code, coords, feature_db)
    assert v.start == coding_start
    assert v.chrom == chrom
    assert v.g.posedit.pos.start.base == genomic_start


# https://towardsdatascience.com/pytest-for-data-scientists-2990319e55e6
signature = 'code, feature_db, is_delta, is_unique'
testdata = [
    ('NM_000546.6:c.(?_560-1)_(672+1_?)del', db, True, True),
    ('NM_015015.3:c.2441+1G>A', db, False, None),
    ('NM_000546.6:g.(123_234567)_(345678_?)del', db, True, None),
        # NM_000546.6:c.(4071+1_4072-1)_(5154+1_5155-1)del
    
        # NM_000546.6:c.(?_560-1)_(672+1_?)del
    
        # g.(?_234567)_(345678_?)del           -- deleted exon is (234567, 345678)
        # c.(4071+1_4072-1)_(5154+1_5155-1)del -- deleted exon is (4072, 5154)
    ]


@pytest.mark.parametrize(signature, testdata)
def test_delta(code, feature_db, is_delta, is_unique):
    ed = ExonDelta(code, feature_db)
    assert ed.is_delta == is_delta
    assert ed.is_unique == is_unique
        


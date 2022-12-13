import json
import pytest

from cdot.hgvs.dataproviders import JSONDataProvider
from gffutils import FeatureDB

from primer4.models import Variant, ExonDelta, Template
from primer4.utils import sync_tx_with_feature_db


@pytest.fixture(scope='session')
def path_to_config(request):
    return request.config.getoption("--config")


@pytest.fixture(scope='session')
def params(path_to_config):
    with open(path_to_config, 'r') as file:
        params = json.load(file)
    return params
'''
"params()" wont work bc/

> Calling a fixture function directly, as opposed to request them in a test function, is deprecated. -- https://docs.pytest.org/en/stable/deprecations.html#calling-fixtures-directly

Related:

https://github.com/pytest-dev/pytest/issues/3950
'''


@pytest.fixture(scope='session')
def hdp(params):
    # assert params['data']['coordinates'] == 1
    # return JSONDataProvider(params['data']['coordinates'])
    return JSONDataProvider([params['data']['coordinates']])


@pytest.fixture(scope='session')
def db(params):
    return FeatureDB(params['data']['annotation'], keep_order=True)


def test_print_name(path_to_config):
    print(f"\ncommand line param (name): {path_to_config}")


def test_presence(params):
    assert 'annotation' in params['data'].keys()


def test_variant(db):
    # v = Variant('NM_000546.6:c.215C>G', hdp, db)
    return db['rna-NM_000546.6']


# Test multiple variants using fixture parametrization:
# - https://towardsdatascience.com/pytest-for-data-scientists-2990319e55e6
# - https://stackoverflow.com/questions/61475337/refreshing-pytest-fixtures-in-first-test-during-custom-scenario-runner
partial_fn_signature = 'version, code, coding_start, chrom, genomic_start'
testdata = [
    ('hg19', 'NM_000546.6:c.215C>G', 215, 'NC_000017.10', 7579472),
    ('hg19', 'NM_000546.6:c.672+3C>G', 672, 'NC_000017.10', 7578174),
    ('hg19', 'NM_015015.3:c.2441+1G>A', 2441, 'NC_000019.9', 5137688)
    ]
# TODO: Add hg38 validation data


# @pytest.mark.parametrize('number1,number2', [(1, 2), (3, 4)])
@pytest.mark.parametrize(partial_fn_signature, testdata)
def test_variant2(version, code, coding_start, chrom, genomic_start, hdp, db, params):
    v = Variant(code, hdp, params['version'])
    tmp = Template(v, db)
    
    if version == params['version']:
        assert v.start == coding_start
        assert v.chrom == chrom
        assert v.g.posedit.pos.start.base == genomic_start
        
        assert tmp.relative_pos(tmp.start) == 0
        assert tmp.relative_pos(tmp.end) == len(tmp)
        assert tmp.invert_relative_pos(tmp.relative_pos(tmp.start)) == tmp.start


# NM_000546.6:c.(4071+1_4072-1)_(5154+1_5155-1)del
# g.(?_234567)_(345678_?)del           -- deleted exon is (234567, 345678)
# c.(4071+1_4072-1)_(5154+1_5155-1)del -- deleted exon is (4072, 5154)
signature = 'code, is_delta, is_unique'
testdata = [
    ('NM_000546.6:c.(?_560-1)_(672+1_?)del', True, True),
    ('NM_015015.3:c.2441+1G>A', False, None),
    ('NM_000546.6:g.(123_234567)_(345678_?)del', True, None),
    ]


@pytest.mark.parametrize(signature, testdata)
def test_delta(code, is_delta, is_unique, db):
    ed = ExonDelta(code, db)
    assert ed.is_delta == is_delta
    assert ed.is_unique == is_unique
      

def test_version_mismatch_tx(db):
    '''
    Test version mismatch between query and feature database
    '''
    assert 'NM_005585.5' == sync_tx_with_feature_db('NM_005585.4', db)



import re
from uuid import uuid4

import click
from gffutils.feature import Feature
import hgvs
from hgvs.assemblymapper import AssemblyMapper
import pyfaidx
from pysam import VariantFile

from primer4.methods import sanger, qpcr, mrna
from primer4.space import pythonic_boundaries
from primer4.utils import (
    convert_chrom,
    gc_map,
    load_variation,
    log,
    manual_c_to_g,
    sync_tx_with_feature_db)


class Variant():
    '''
    hdp = JSONDataProvider(['cdot-0.2.1.refseq.grch37_grch38.json.gz'])
    v = Variant('NM_000546.6:c.215C>G', hdp)

    HGVS syntax examples:

    - NM_000546.6:c.215C>G
    - NM_000546.6:c.215_250del
    - NM_000546.6:c.215_250del

    - g.(?_234567)_(345678_?)del           -- deleted exon is (234567, 345678)
    - c.(4071+1_4072-1)_(5154+1_5155-1)del -- deleted exon is (4072, 5154)
    - (6278_6438+69)_(7310-43_7575)del)

    - NM_015015.3:c.2441+1G>A
    '''
    def __init__(self, code, coords, feature_db):
        self.code = code
        self.data = self.parse_code(code)
        self.g = None
        self.chrom = None
        self.start = self.data.posedit.pos.start.base
        self.start_offset = self.data.posedit.pos.start.offset
        self.end = self.data.posedit.pos.end.base
        self.end_offset = self.data.posedit.pos.end.offset        
        self.tx = self.data.ac

        if self.start_offset:
            self.is_coding = False
        else:
            self.is_coding = True

        self.g = self.map_to_genomic(coords)
        self.chrom = self.g.ac
        self.g_start = self.g.posedit.pos.start.base
        self.g_end = self.g.posedit.pos.end.base

        # Overwrite tx version to sync w/ annotation; v.tx and v.data.ac
        self.tx = sync_tx_with_feature_db(self.tx, feature_db)

    def __repr__(self):
        return self.code

    def parse_code(self, code):
        hp = hgvs.parser.Parser()
        try:
            return hp.parse_hgvs_variant(code)
        except hgvs.exceptions.HGVSParseError:
            if ed := is_exon_deletion(code):
                return ed
            else:
                # click.echo(log('Invalid HGVS syntax, exit.'))
                print('Variant cannot be HGVS-parsed!')
                return None

    def map_to_genomic(self, tx_map, genome_version='GRCh37', method='splign'):
        am = AssemblyMapper(
            tx_map,
            assembly_name=genome_version,
            alt_aln_method=method,
            replace_reference=True)
        return am.c_to_g(self.data)


class ExonDelta():
    def __init__(self, code, feature_db):
        self.code = code
        self.tx, self.is_delta, self.data = self.is_exon_delta(code, feature_db)
        
        if not self.data:
            self.is_unique = None
        else:
            self.is_unique = True if len(self.data) == 1 else False
        
        if self.is_unique:
            x = list(self.data)[0]
            self.chrom = x.chrom
            self.g_start = x.start
            self.g_end = x.end


    def __repr__(self):
        return self.data.__repr__()

    def is_exon_delta(self, code, feature_db):
        '''
        NM_000546.6:g.(?_234567)_(345678_?)del
        NM_000546.6:g.(123_234567)_(345678_?)del
        NM_000546.6:c.(4071+1_4072-1)_(5154+1_5155-1)del
    
        NM_000546.6:c.(?_560-1)_(672+1_?)del
    
        g.(?_234567)_(345678_?)del           -- deleted exon is (234567, 345678)
        c.(4071+1_4072-1)_(5154+1_5155-1)del -- deleted exon is (4072, 5154)
    
        is_exon_deletion(code, db)
        # <Feature exon (NC_000017.10:7578177-7578289[-]) at 0x7f9c32908d60>
        '''
        num = r'[\?]?\d*?[+-]?\d*?'
        cap = r'(\d*)[+-]?\d*?'  # cap .. capture
        pattern = rf'(.*):([gc])\.\({num}_{cap}\)_\({cap}_{num}\)(.*)'
    
        m = re.match(pattern, code)
        
        try:
            # Try to parse code
            tx = m.group(1)
            tx = sync_tx_with_feature_db(tx, feature_db)
        except AttributeError:
            # Not an exon delta;
            # AttributeError: 'NoneType' object has no attribute 'group'
            tx = code.split(':')[0]
            tx = sync_tx_with_feature_db(tx, feature_db)
            return tx, False, set()

        is_coding = True if m.group(2) == 'c' else False
        start     = int(m.group(3))  # coding start coord
        end       = int(m.group(4))  # ... end ...

        if is_coding:
            start = manual_c_to_g(tx, start, feature_db)
            end = manual_c_to_g(tx, end, feature_db)
            if start > end: start, end = end, start
            # Otherwise gffutils does not find a .region()

        A = set(feature_db.children(f'rna-{tx}', featuretype='exon', order_by='start'))
        
        B = set(feature_db.region(region=(list(A)[0].chrom, start, end), featuretype='exon'))
        AB = A.intersection(B)
        
        return tx, True, AB


class ExonSpread():
    def __init__(self, transcript, exon1, exon2):
        self.tx = transcript
        self.exon1 = exon1
        self.exon2 = exon2
        self.data = (transcript, exon1, exon2)

    def __repr__(self):
        return f'{self.tx}, (exon {self.exon1} and {self.exon2})'


class SingleExon():
    def __init__(self, transcript, exon):
        self.tx = transcript
        self.exon = exon
        self.data = (transcript, exon)

    def __repr__(self):
        return f'{self.tx}, exon {self.exon}'

'''
Template:

- gets tx and retrieves/ masks variants
- gets either a variant or and exon (chrom, start, end)
- applies sanger to var or qpcr/ splice to exon

(decouple input from fn applied, bc/ n:n mapping)

    t = Template(v, db)
    t.relative_pos(t.start)
    # 0

    Needs to interact w/

    - variant
    - manual parse

    Needs to contain:

    - tx coords and relative translation
    - target position in tx
    - padded sequences and other masks (cannot think of any but there might be)
    - SNVs

    Then

    t = Template(...)
    constraints = t.search_for('sanger')
    design(t.masked_sequence, constraints, params['sanger'])

    Template(v, db)
    Template(ed, db)
'''
class Template():
    def __init__(self, mutation, feature_db, featuretype='exon'):
        self.data = mutation
        self.type = type(mutation)  # Template(v, db).type == Variant
        self.feat = feature_db[f'rna-{mutation.tx}']
        self.region = list(feature_db.region(
            region=(self.feat.chrom, self.feat.start, self.feat.end),
            featuretype=featuretype))
        self.start, self.end = pythonic_boundaries(self.feat)
        self.mask = set()
        self.methods = {
            'sanger': sanger,
            'qpcr':   qpcr,
            'mrna':   mrna,
        }
        self.accepted = {
            'sanger': Variant,
            'qpcr': ExonDelta,
            'mrna': ExonDelta,
        }
        self.variation = {}
        self.c_to_g, self.g_to_c = gc_map(mutation.tx, feature_db)
        # TODO: Maybe duplicate code here to self.feat
        self.mrna = ()

    def __repr__(self):
        return self.feat.__repr__()

    def __len__(self):
        return len(self.feat)

    # TODO: apply mask as fn
    def get_sequence(self, genome):
        s = genome[self.feat.chrom][self.start:self.end].__str__()
        # assert len(s) == len(self.feat)
        return s

    def relative_pos(self, n):
        # Primer3 needs positions relative to sequence (when masking etc.)
        # Turns genomic coordinate into relative one
        return n - self.start
    
    def apply(self, fn, feature_db, params, *args, **kwargs):
        # Check that we apply the right fn to the right data type
        # TODO:
        # assert self.type == self.accepted.get(fn)
        # Apply fn
        return self.methods[fn](self, feature_db, params, *args, **kwargs)

    def load_variation_(self, databases):
        self.mask = load_variation(self.feat, databases)
        return None

    # def mask_sequence(self, genome, mask='N', unmasked=''):
    #     s = self.get_sequence(genome)
        
    #     if unmasked:
    #         masked = ''.join([mask if ix in self.mask else unmasked for ix, i in enumerate(s)])
    #     else:
    #         masked = ''.join([mask if ix in self.mask else i for ix, i in enumerate(s)])
    #     return masked.upper()


class PrimerPair():
    '''
    https://stackoverflow.com/questions/1305532/convert-nested-python-dict-to-object/9413295#9413295
    '''
    def __init__(self, d):
        self.name = uuid4().__str__()
        self.data = d

        for a, b in d.items():
            if isinstance(b, (list, tuple)):
               setattr(self, a, [PrimerPair(x) if isinstance(x, dict) else x for x in b])
            else:
               setattr(self, a, PrimerPair(b) if isinstance(b, dict) else b)

    def __repr__(self):
        return f'{self.fwd.start}-{self.fwd.end}:{self.rev.start}-{self.rev.end}, loss: {self.penalty}'

    def to_g(self, template):
        '''
        Translate primer coords into genomics ones.
        '''
        fwd_start = template.feat.start + self.fwd.start
        fwd_end   = template.feat.start + self.fwd.end
        rev_start = template.feat.start + self.rev.start
        rev_end   = template.feat.start + self.rev.end    
    
        return [fwd_start, fwd_end, rev_start, rev_end]

    def to_c(self, template):
        '''
        Translate primer coords into coding ones.

        TODO: Can take genomic position and query template for the start of
        the closest exon, or pass spanned exon and calc distance then write
        coding_start - difference.

        If its on the exon, have this position, else return the start of the
        exon?
        '''
        # TODO: Does not work if non-coding position
        return [template.g_to_c[i] for i in self.to_g(template)]
        '''
        TODO: Get the closest coding position, then do +- x.
        '''
        # qry = 7579200
        # nearest_g = min(tmp.g_to_c.keys(), key=lambda x: abs(x - qry))
        # nearest_c = tmp.g_to_c[nearest_g]
        # str(qry - nearest_g)  # -112
        # TODO: start or end of primer reference here ie -132 oder -112

    def save(self, fp):
        with open(fp, 'w+') as out:
            for i in ['fwd', 'rev']:
                out.write(f'>{self.name}.{i}\n{self.data[i]["sequence"]}\n')

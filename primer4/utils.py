import datetime
from difflib import get_close_matches
from itertools import chain, zip_longest
from math import floor
from pathlib import Path
import pdb
# pdb.set_trace()
import sys

import click
from hgvs.parser import Parser
# https://github.com/biocommons/hgvs
import primer3
from pysam import VariantFile


# https://stackoverflow.com/questions/1270951/how-to-refer-to-relative-paths-of-resources-when-working-with-a-code-repository
fn = Path(__file__).parents[1] / 'data/chrom_names.csv'
chrom_names = {}
with open(fn, 'r') as file:
    for line in file:
        k, v = line.strip().split(',')
        chrom_names[k] = v


def convert_chrom(chrom):
    return chrom_names.get(chrom)


def infer_coordinates(variant, db):
    '''
    Ah!
    
    https://www.biostars.org/p/104870/
    
    - Exons = gene - introns
    - CDS   = gene - introns - UTRs
    
    and by arithmetic:
    
    - CDS   = Exons - UTRs
    
    So need to look at CDS, not exons.
    
    gffutils assigns IDs to features; for CDS it calls the first CDS in a transcript like "cds-NP_612468.1" (note the prefix) and the next cds-NP_612468.1_1 and so on. Tested for:
    
    - NM_000546.6
    - NM_138459.5
    '''
    hp = Parser()
    try:
        v = hp.parse_hgvs_variant(variant)
    except hgvs.exceptions.HGVSParseError:
        click.echo(log('Invalid HGVS syntax, exit.'))
        sys.exit(-1)

    # NM_000350.3:c.4773+3A>G
    # v.posedit.pos.start
    # BaseOffsetPosition(base=4773, offset=3, ...)
    c_pos = v.posedit.pos.start.base
    # coding pos, base (e.g. 4773 from 4773+3)
    offset = v.posedit.pos.start.offset
    if offset:
        click.echo(log('Non-coding variant'))
    name = v.ac

    # Generate a list of valid IDs so we can validate input:
    IDs = set(i.id.replace('rna-', '') for i in db.features_of_type('mRNA'))

    if not name in IDs:
        m = get_close_matches(name, IDs)
        q = name.split('.')[0]
        
        ID_version = {k: v for k, v in [i.split('.') for i in IDs]}
        v = ID_version[q]
        name = f'{q}.{v}'
        click.echo('\n' + log('Warning!\n'))
        click.echo(f'Transcript ID not found, version mismatch? Will use {name} instead.')
        click.echo(f'Not what you want? How about: {m[0]}, {m[1]} or {m[2]}?')

    chromosome = db[f'rna-{name}'].chrom


    coding = []
    for i in db.children(f'rna-{name}', featuretype='CDS', order_by='start'):
        coding.append(i)
    
    strand = db[f'rna-{name}'].strand
    l = 0
    if strand == '-':
        coding = coding[::-1]  # reverse CDS order
        
        for i in coding:
            if l + len(i) >= c_pos:
                break
            else:
                l += len(i)
        # i.end holds genomic coordinate
        g_pos = i.end - (c_pos - l) + 1
    
    else:
        for i in coding:
            if l + len(i) >= c_pos:
                break
            else:
                l += len(i)
        g_pos = i.start + c_pos - 1


    if offset:
        # The offset is in relation to the coding sequence, eg "+3" means
        # three bases downstream of the previous coding sequence.
        if strand == '-':
            g_pos -= offset
        else:
            g_pos += offset

    click.echo(log(f'Variant on chromosome {chromosome}, g. position {g_pos}'))
    return name, chromosome, g_pos


def retrieve_exon_around_variant(db, name, chromosome, g_pos):
    '''
    Return exon boundaries either for a specified (transcript, exon) pair or
    a genomic position.

    Returns start and end coords in "Python space" so slicing works as expected.

    !grep 'exon-NM_000546.6-4' GRCh37_latest_genomic.gff
    NC_000017.10    BestRefSeq  exon    7579312 7579590 .   -   .   ID=exon-NM_000546.6-4;Parent=rna-NM_000546.6;Dbxref=GeneID:7157,Genbank:NM_000546.6,HGNC:HGNC:11998,MIM:191170;gbkey=mRNA;gene=TP53;product=tumor protein p53%2C transcript variant 1;tag=RefSeq Select;transcript_id=NM_000546.6

    import gffutils
    db = gffutils.FeatureDB('hg19-p13_annotation.db', keep_order=True)
    ex = db['exon-NM_000546.6-4']
    ex
    # <Feature exon (NC_000017.10:7579312-7579590[-]) at 0x7fb001c14b20>
    '''
    # Get all exons that cover the genomic position and select the one
    # corresponding to the transcript <name>
    g = db.region(region=(chromosome, g_pos, g_pos + 1), featuretype='exon')
    exons = [i for i in g]

    for ex in exons:
        # An exon can be part of multiple transcripts; choose the exon
        # annotation particular to the transcript of interest. We assume
        # that the exons of any particular transcript do not overlap.
        if name in ex.id:
            return ex

    # A variant can be intronic, in which case we don't find any exon
    return None


def retrieve_exon(db, name, exon):
    return db[f'exon-{name}-{exon}']


def pythonic_exon_boundaries(exon):
    # Get Python coords for intuitive sclicing later
    if exon.strand == '-':
        start = exon.start - 1
        end = exon.end
    else:
        start = exon.start
        end = exon.end + 1
    return start, end


def exon_context(name, exon, db, params):
    # left in right out OR left out right in OR both in
    mn, mx = params['size_range_qPCR']
    binding = params['binding_site']

    ex = retrieve_exon(db, name, exon)
    start, end = pythonic_exon_boundaries(ex)

    # TODO
    if mn > len(ex):
        pass
    else:
        pass

    left  = start - binding
    right =   end + binding
    s = genome[ex.chrom][left:right].__str__()
    
    pass



def variant_context(name, genome, chromosome, g_pos, db, params):
    '''
    Cases:

    [x] exon spanned
    [x] exon not spanned
    [ ] exon lost and qPCR
    '''
    # Primers need 18-30 nt and then we leave another 30 for Sanger "burn-in";
    # so we'll introduce a padding parameter
    # https://eu.idtdna.com/pages/support/faqs/what-is-the-optimal-length-of-a-primer-
    pad = params['burnin_sanger']
    binding = params['binding_site']  # minimum search space for each primer
    # Minimum, maximum amplicon size
    _, mx = params['size_range_PCR']

    if ex := retrieve_exon_around_variant(db, name, chromosome, g_pos):
        required = binding + pad + len(ex) + pad + binding
    else:
        # No exon annotated for the variant
        required = 1e10
        # Some insanely large value to trigger variant-centered primer design

    click.echo(log('Extracting target sequence; depending on chromosome size might take a while'))
    if mx > required:  # as in the sequence
        click.echo(log('Will span exon'))
        # Mask the entire exon and try to find primers spanning the exon
        # Then fill up left and right to maximum amplicon size
        fill = floor((mx - required) / 2)
        
        start, end = pythonic_exon_boundaries(ex)
        left  = start - pad - binding - fill
        right =   end + pad + binding + fill

        '''
        Template:

        left (fill, binding) (pad) start (exon) end (pad) (fill, binding) right
        '''
        s = genome[chromosome][left:right].__str__()
        
        s1 = s[:binding + fill].upper()
        s2 = 'N' * (len(ex) + (2 * pad))
        s3 = s[len(s1 + s2):].upper()
        template = s1 + s2 + s3
        assert len(s) == len(template)

        # Where is it acceptable to search for left and right primers, resp.?
        # SEQUENCE_PRIMER_PAIR_OK_REGION_LIST
        # https://primer3.ut.ee/primer3web_help.htm
        r = fill + binding
        left_search  = (1, r)
        right_search = (len(template) - r, r)

    else:
        click.echo(log('Could not span exon'))
        # Now put an amplicon-sized window across the mutation, such that the
        # mutation has enough space to both sides. Primer3 can later find 
        # primers in this window subject to the (mn, mx) constraints
        fill = floor((mx - (2 * pad)) / 2)
    
        left  = g_pos - pad - fill
        right = g_pos + pad + fill
    
        s = genome[chromosome][left:right].__str__()
    
        s1 = s[:fill].upper()
        s2 = 'N' * (2 * pad + 1)  # add + 1 for g_pos
        s3 = s[len(s1 + s2):].upper()
        template = s1 + s2 + s3
        assert len(s) == len(template)

        r = fill
        left_search  = (1, r)
        right_search = (len(template) - r, r)

    return template, (left, right), (left_search, right_search)


def common_variants(template, chromosome, boundaries, variants):
    left, right = boundaries
    # Don't put primers in the vicinity of the mutation
    # common, mask = [], []
    mask = set()
    # Don't put them across positions with common variants
    for i in variants.fetch(chromosome, left, right):
         if i.info.get('COMMON'):
            # print(i)
            # common.append(i.ref)
            mask.add(i.pos - left - 1)
            # print(i.ref, s[i.pos - left - 1])
            # print(genome[chromosome][i.pos - 1], i.ref)
            # Not sure why, but we have to subtract 1 to obtain the same letter
    
    masked = ''.join(
        ['N' if ix in mask else i for ix, i in enumerate(template)])
    
    click.echo(log('Template ("x" means no primers here):\n'))
    print(''.join(['x' if i == 'N' else '-' for i in masked]))
    return masked


def design_primers(method, params, masked, constraints):
    '''
    # 'SEQUENCE_EXCLUDED_REGION'
    # https://primer3.org/manual.html#SEQUENCE_EXCLUDED_REGION
    
    # 'SEQUENCE_INCLUDED_REGION'
    
    # 'PRIMER_PRODUCT_SIZE_RANGE'
    # https://github.com/libnano/primer3-py/issues/18
    '''
    size_range = params[f'size_range_{method}']

    # https://primer3.ut.ee/primer3web_help.htm
    # SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=100,50,300,50 ; 900,60,, ; ,,930,100
    (left_start, left_len), (right_start, right_len) = constraints
    only_here = list(chain(*constraints))
    # print(only_here)
    # pdb.set_trace()

    spec =  [
        {
            'SEQUENCE_TEMPLATE': masked,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': only_here,
        },
        {
            'PRIMER_NUM_RETURN': params['n_return'],
            'PRIMER_MIN_SIZE': params['size_min'],
            'PRIMER_OPT_SIZE': params['size_opt'],
            'PRIMER_MAX_SIZE': params['size_max'],
    
            'PRIMER_MIN_TM': params['tm_min'],
            'PRIMER_OPT_TM': params['tm_opt'],
            'PRIMER_MAX_TM': params['tm_max'],
            'PRIMER_MIN_GC': params['GC_min'],
            'PRIMER_MAX_GC': params['GC_max'],
    
            'PRIMER_MAX_POLY_X': params['homopolymer_max_len'],
            'PRIMER_MAX_END_GC': params['3prime_max_GC'],
    
            'PRIMER_MAX_NS_ACCEPTED': params['Ns_max'],
            
            'PRIMER_PRODUCT_SIZE_RANGE': [size_range],
    
            # defaults, here to be explicit
            'PRIMER_SALT_MONOVALENT': params['salt_monovalent'],
            'PRIMER_SALT_DIVALENT': params['salt_divalent'],
            'PRIMER_DNTP_CONC': params['conc_dNTP'],
            'PRIMER_DNA_CONC': params['conc_DNA'],
        }
    ]

    # https://libnano.github.io/primer3-py/quickstart.html#workflow
    design = primer3.bindings.designPrimers(*spec)
    return design


def log(message):
    now = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
    return f'[{now}]\t{message}'


def parse_design(design, n, pname, g_pos, chromosome, masked):
    
    primers = {}
    for i in range(n):
        
        fwd_start, fwd_len = design[f'PRIMER_LEFT_{i}']
        rev_start, rev_len = design[f'PRIMER_RIGHT_{i}']
        
        # c .. candidate
        primers[f'{pname}::c{i}'] = {
            'fwd': {
                'start': fwd_start,
                'end': fwd_start + fwd_len,
                # 'sanity': template[fwd_start:fwd_start + fwd_len],
                'sequence': design[f'PRIMER_LEFT_{i}_SEQUENCE'],
                'Tm': round(design[f'PRIMER_LEFT_{i}_TM'], 2),
            },
            'rev': {
                # That Python 0-based, end-exclusive indexing thing ...
                'start': rev_start - rev_len + 1,
                'end': rev_start + 1,
                # 'sanity': rc(template[rev_start - rev_len + 1:rev_start + 1]),
                'sequence': design[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                'Tm': round(design[f'PRIMER_RIGHT_{i}_TM'], 2),
            },
            'insert': design[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
            'template': masked,
            'penalty': round(design[f'PRIMER_PAIR_{i}_PENALTY'], 4),
            'chromosome': chromosome,
            'genomic_position': g_pos,
        }

    return primers


def manual_c_to_g(tx, c, feature_db):
    '''
    Turn a coding variant into a genomic one

    # NM_000546.6:c.215C>G
    manual_c_to_g('NM_000546.6', 215, db)
    # 7579472
    '''
    coding = []
    for i in feature_db.children(
        f'rna-{tx}', featuretype='CDS', order_by='start'):
        coding.append(i)
    
    strand = feature_db[f'rna-{tx}'].strand
    l = 0
    if strand == '-':
        coding = coding[::-1]  # reverse CDS order
        
        for i in coding:
            if l + len(i) >= c:
                break
            else:
                l += len(i)
        # i.end holds genomic coordinate
        g = i.end - (c - l) + 1
    
    else:
        for i in coding:
            if l + len(i) >= c:
                break
            else:
                l += len(i)
        g = i.start + c - 1
    
    return g


def sync_tx_with_feature_db(tx, feature_db):
    # Generate a list of valid IDs so we can validate input:
    IDs = set(i.id.replace('rna-', '') for i in feature_db.features_of_type('mRNA'))

    if not tx in IDs:
        # m = get_close_matches(name, IDs)
        q = tx.split('.')[0]
    
        ID_version = {k: v for k, v in [i.split('.') for i in IDs]}
        v = ID_version[q]
        name = f'{q}.{v}'
        click.echo('\n' + log('Warning!\n'))
        click.echo(f'Transcript ID not found, version mismatch? Will use {name} instead.')
        # click.echo(f'Not what you want? How about: {m[0]}, {m[1]} or {m[2]}?')
        chromosome = feature_db[f'rna-{name}'].chrom
        return tx

    else:
        return tx


def gc_map2(tx, feature_db):
    '''
    Map coding coordinates to genomic ones for easy lookup.
    '''
    cnt, cum = 0, 0
    map_ = {}
    for i in feature_db.children(
        f'rna-{tx}', featuretype='CDS', order_by='strand'):
        cum += len(i)

        r = range(i.start, i.end+1)
        if i.strand == '-':
            r = reversed(r)
        
        for j in r:
            cnt += 1
            map_[cnt] = j

        assert max(map_.keys()) == cum
    return map_, dict((v, k) for k, v in map_.items())



def gc_map(tx, feature_db):
    '''
    Map coding coordinates to genomic ones for easy lookup.
    '''
    cnt, cum = 0, 0
    map_ = {}

    tx = f'rna-{tx}'
    
    strand = feature_db[tx].strand
    l = list(feature_db.children(tx, featuretype='CDS', order_by='start'))

    if strand == '-':
        l = l[::-1]

    for i in l:  # r .. region
        cum += len(i)
        
        r = range(i.start, i.end+1)
        if i.strand == '-':
            r = reversed(r)

        for j in r:
            cnt += 1
            map_[cnt] = j

    assert max(map_.keys()) == cum
    return map_, dict((v, k) for k, v in map_.items())


def twolists(l1, l2):
    '''
    https://stackoverflow.com/questions/48199961/   how-to-interleave-two-lists-of-different-length
    '''
    return [x for x in chain(*zip_longest(l1, l2)) if x is not None]


def load_variation(feat, databases):
    '''
    feat .. gffutils feature type
    '''
    mask = set()
    for name, db in databases.items():
        variants = VariantFile(db)

        # dbSNP names chromosomes like "NC_000007.13", others like "7"
        if name != 'dbSNP':
            vv = variants.fetch(convert_chrom(feat.chrom), feat.start, feat.end)
        else:
            vv = variants.fetch(feat.chrom, feat.start, feat.end)

        for i in vv:
            # .info.get(...) raises ValueError: Invalid header if not there
            info = dict(i.info)
            
            if name == 'dbSNP':
                if info.get('COMMON'):
                    mask.add(i.pos - feat.start - 1)
            
            elif name == '1000Genomes':
                if info['AF'][0] >= 0.01:
                    mask.add(i.pos - feat.start - 1)
            
            elif name == 'ESP':
                if float(info['MAF'][0]) >= 1:
                    mask.add(i.pos - feat.start - 1)

            else:
                print(f'"{name}" is not a valid variant database')
    return mask


def mask_sequence(seq, var, mask='N', unmasked=''):
    
    if unmasked:
        masked = ''.join([mask if ix in var else unmasked for ix, i in enumerate(seq)])
    else:
        masked = ''.join([mask if ix in var else i for ix, i in enumerate(seq)])
    return masked.upper()


def reconstruct_mrna(tx, feature_db, genome, vardbs):
    exons = {}
    for e in feature_db.children(tx, featuretype='exon', order_by='start'):
        exons[int(e.id.split('-')[-1])] = e

    reconstruction = ''
    coords = []
    segmentation = []

    # print(exons)
    for k in sorted(exons.keys()):
        ex = exons[k]
        
        # seq = ex.sequence(genome).upper()  # accounts for strand

        #seq = sequence[ex.start:ex.end+1]
        #print(sequence[:10], ex.start, ex.end)
        
        seq = ex.sequence(genome).upper()
        if vardbs:
            var = load_variation(ex, vardbs)
            seq = mask_sequence(seq, var) 

        reconstruction += seq
        # Validated manually using screed that this considers strand, ie for "-"
        # we get the revcomp sequence.
    
        pos = list(range(ex.start, ex.end+1))
        if ex.strand == '-':
            pos = list(reversed(pos))
            # +1 bec intervals INCLUDE the last position, eg 7573008 below, but
            # the reversed fn() excludes it:
            # <Feature exon (NC_000017.10:7571739-7573008[-]) at 0x7fb7fe96eca0>
            # list(reversed(range(ex.start, ex.end)))[0]   is 7573007
            # list(reversed(range(ex.start, ex.end+1)))[0] is 7573008 
        # print(len(pos), len(ex), len(seq))
        assert len(pos) == len(ex) == len(seq)
    
        segmentation.extend([k] * len(seq))
        coords.extend(pos)
    
    return reconstruction, exons, coords, segmentation


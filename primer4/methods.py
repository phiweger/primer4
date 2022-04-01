from primer4.space import context
from primer4.utils import twolists#, reconstruct_mrna


def sanger(template, feature_db, params): 
    '''
    mask_sanger(v, t, params)
    '''
    mn, mx = params['size_range_PCR']
    burnin = params['burnin_sanger']

    pad = burnin * 2 + 100
    # x2 bc/ burnin on both sides of variant, 100 is space to search primers in 
    variant = template.data
    ex = context(variant, feature_db, 'exon')
    
    if ex and (len(ex) + pad < mx):
        # We can span the exon
        # exon boundary + pad, rest primer
        lb = ex.start - burnin  # left boundary
        rb = ex.end + burnin    # right boundary
    else:
        # Cannot span exon
        # center + pad, rest primer
        
        lb = variant.g_start - burnin
        rb = variant.g_end + burnin

    # Primer3 takes (start, length) constraints
    rlb = template.relative_pos(lb)
    rrb = template.relative_pos(rb)
    return {
        'only_here': ((0, rlb), (rrb, len(template) - rrb)),
        'size_range': (mn, mx)
    }


def qpcr(template, feature_db, params):
    '''
    qpcr(tmp, db, params, 5)
    '''
    mn, mx = params['size_range_qPCR']

    exons = list(feature_db.children(
        template.feat.id, featuretype='exon', order_by='start'))
    introns = list(feature_db.interfeatures(exons))
    l = twolists(exons, introns)
    ix = [int(i.id.split('-')[-1]) if i.id else None for i in l]
    mid = ix.index(template.data.exon)
    left = mid - 1 
    right = mid + 1

    # left
    rlb = template.relative_pos(l[left].start)  # rlb .. relative left bound
    rmb = template.relative_pos(l[mid].start)  # rmb .. mid
    cl = {
        'only_here': (
            (rlb, len(l[left])),
            (rmb, len(l[mid]))
            ),
        'size_range': (mn, mx)
    }
    # right

    rrb = template.relative_pos(l[right].start)  # rrb .. relative right bound
    cr = {
        'only_here': (
            (rmb, len(l[mid])),
            (rrb, len(l[right]))
            ),
        'size_range': (mn, mx)
        }
    return cl, cr


def mrna(template, feature_db, params):
    if not template.mrna:
        raise ValueError('Please reconstruct the mRNA first')
    
    mrna, exons, coords, labels = template.mrna

    here = []
    for exon in template.data.data[1:]:
        x = [i for i, j in enumerate(labels) if j == exon]
        mn, mx = min(x), max(x)
        here.append((mn, mx-mn))

    constraints = {
        'only_here': tuple(here),
        'size_range': tuple(params['size_range_mRNA'])
    }

    return constraints


'''
def mrna(template, feature_db, params):
    mrna, exons, coords, labels = reconstruct_mrna(template.feat.id, feature_db)

    here = []
    for exon in template.data.data[1:]:
        x = [i for i, j in enumerate(labels) if j == exon]
        mn, mx = min(x), max(x)
        here.append((mn, mx-mn))

    constraints = {
        'only_here': tuple(here),
        'size_range': tuple(params['size_range_qPCR'])
    }

    return constraints
'''

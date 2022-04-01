def context(variant, feature_db, featuretype='exon'):
    '''
    context(v, db, 'exon')
    '''
    chrom = variant.chrom
    start = variant.g.posedit.pos.start.base
    end = variant.g.posedit.pos.end.base
    
    r = feature_db.region(region=(chrom, start, end), featuretype=featuretype)

    for feat in r:
        # An exon can be part of multiple transcripts; choose the exon
        # annotation particular to the transcript of interest. We assume
        # that the exons of any particular transcript do not overlap.
        if variant.tx in feat.id:
            return feat
    
    # A variant can be intronic, in which case we don't find any exon
    return None


def neighbor(exon, db, direction=-1):
    '''
    v = Variant('NM_000546.6:c.215C>G', hdp, db)
    ex = context(v, db, 'exon')
    ex.id
    # 'exon-NM_000546.6-4'
    neighbor(ex, db, -1).id
    # 'exon-NM_000546.6-3'

    The strand is considered.
    '''
    _, tx, n = exon.id.split('-')  # exon-NM_000546.6-4
    return db[f'exon-{tx}-{int(n) + direction}']


def pythonic_boundaries(feature):
    # Get Python coords for intuitive sclicing later
    if feature.strand == '-':
        start = feature.start - 1
        end = feature.end
    else:
        start = feature.start
        end = feature.end + 1
    return start, end

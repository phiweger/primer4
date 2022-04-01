import gffutils

from primer4.utils import infer_coordinates


'''
# Genomic coords
fp_rna = 'data/GRCh37_latest_rna.fna.gz'
rna = Fasta(fp_rna)
tx = db['rna-NM_000546.6']
a = tx.sequence(fp)[:10]
# <Feature mRNA (NC_000017.10:7571739-7590808[-]) at 0x7ff3806406d0>
# Note: "-" strand
b = rc(genome[chromosome][tx.start:tx.end].__str__())[:10]
assert a == b

ex = db['exon-NM_000546.6-11']
a = ex.sequence(fp)[:10]
b = rc(genome[chromosome][ex.start:ex.end].__str__())[-10:]
assert a == b
'''


'''

# TODO: Add as fixture
db = gffutils.FeatureDB(fp_annotation, keep_order=True)


genome = Fasta(fp_genome)
# Variant is on chromosome 17
assert genome['NC_000017.10'][87:100].__str__() == 'CCACGACCAACTC'


# 0- or 1-based coordinate systems 
# https://www.biostars.org/p/84686/
# https://plastid.readthedocs.io/en/latest/concepts/coordinates.html
# So, on if the feature is on the "-" strand, start - 1, else end + 1, bc/
# gff is end-inclusive while python is excludes the end.
a = rc(genome[chromosome][ex.start-1:ex.end].__str__())
b = ex.sequence('data/GRCh37_latest_genomic.fna.gz', use_strand=True)
assert a == b


def test_infer_coordinates():
    name, chromosome, g_pos, c_pos = infer_coordinates('NM_138459.3', db)
    assert name == 'NM_138459.6'
    assert chromosome == 'NC_000017.10'
    assert g_pos == 7579472
    assert c_pos == 215
'''
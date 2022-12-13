import re

import pytest

from primer4.design import parse_blast_btop


# https://www.ncbi.nlm.nih.gov/books/NBK569862/
signature = 'query, split, aln'
testdata = [
    ('3G-13', [3, 'G-', 13], '...|.............'),
    ('7AG10', [7, 'AG', 10], '.......|..........'),
    ('7A-10', [7, 'A-', 10], '.......|..........'),
    ('6-G-A10', [6, '-G', '-A', 10], '......||..........'),
    ('6G--A10', [6, 'G-', '-A', 10], '......||..........'),
    ('6G-A-10', [6, 'G-', 'A-', 10], '......||..........'),
    ('6-GA-10', [6, '-G', 'A-', 10], '......||..........'),
    ]

@pytest.mark.parametrize(signature, testdata)
def test_parse_blast_btop(query, split, aln):
     a, b = parse_blast_btop(query, debug=True)
     assert a == aln and b == split

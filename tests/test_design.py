import re

from primer4.design import parse_blast_btop


def test_parse_blast_btop():
    s = '3GA--13'
    assert parse_blast_btop(s) == '...|||.............'

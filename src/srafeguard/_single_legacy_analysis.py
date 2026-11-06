# Verbatim scientific helpers extracted from validated WGS_CONTA_2026 analysis.py.

import re

CIGAR_RE=re.compile(r'(\d+)([MIDNSHP=X])')

def cigar_stats(cigar):
    q=a=r=0
    for n,op in CIGAR_RE.findall(cigar):
        n=int(n)
        if op in 'MI=X': q+=n
        if op in 'MID=X': a+=n
        if op in 'MDN=X': r+=n
    return q,a,r

def parse_nm(f):
    for x in f[11:]:
        if x.startswith('NM:i:'): return int(x.split(':')[2])
    return 0

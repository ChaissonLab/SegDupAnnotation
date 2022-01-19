#!/usr/bin/env python
import pysam
import sys

def GetCDSFromName(name):
    vals=name.split("|")
    for v in vals:
        if "CDS:" in v:
            pair=[int(i)-1 for i in v[4:].split("-")]
            return pair
    return None

af=pysam.AlignmentFile(sys.argv[1])

for rec in af.fetch():
    name=rec.query_name
    cds = GetCDSFromName(name)
    pairs=rec.get_aligned_pairs()
    minStartDiff = 1000000000
    minEndDiff   = 1000000000
    if cds is None:
        continue
    start=cds[0]
    end=cds[1]
    si=None
    ei=None
#    import pdb
#    pdb.set_trace()
    for i in range(0,len(pairs)):
        if pairs[i][0] is None:
            continue
        if abs(pairs[i][0]-start) < minStartDiff:
            si=i
            minStartDiff = abs(pairs[i][0]-start)
        if abs(pairs[i][0]-end) < minEndDiff:
            ei=i
            minEndDiff = abs(pairs[i][0]-end)
            if minEndDiff == 0:
                import pdb
                pdb.set_trace()
                print(name)
    if si is not None and ei is not None:
        print("{}\t{}\t{}".format(minStartDiff, minEndDiff, name))
        


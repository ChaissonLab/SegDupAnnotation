import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-lrt",help="lrt threshold",default=2)
parser.add_argument("-r",help="collapsed region")

args = parser.parse_args()


ratio=[]
for line in sys.stdin:
    line=line.rstrip()
    vals=line.split()
    ref=vals[8]
    counts = [(int(vals[2]), 'A'), (int(vals[3]),'C'), (int(vals[4]),'G'), (int(vals[5]),'T')]
    counts.sort()
    if counts[-1][0]==0 or counts[-2][0]==0:
        continue
    ratio.append(counts[-1][0]/(counts[-2][0]+counts[-1][0]))
    #if abs(ratio-1)<0.4:
     #   print(line)
#    print("\t".join([line,str(ratio)]))
r=np.array(ratio)
#print(len(r[np.where(r>0.609)] )) 
lrt=(len(r[np.where(r>0.609)]) - len(r[np.where(r>0.711)] )) / max(1, (len(r[np.where(r<0.551)  ])  - len(r[np.where(r<0.449)] )))

if lrt>int(args.lrt):
    print(args.r+"\tpass\t"+str(lrt))
else:
    print(args.r+"\tfail\t"+str(lrt))

#!/usr/bin/env python
import sys
bedFile = open(sys.argv[1])
covFile = open(sys.argv[2])
meanCov = float(sys.argv[3])

cov = {}
for line in covFile:
    vals=line.split()
    if vals[0] not in cov:
        cov[vals[0]] = []
        sys.stderr.write("reading " + vals[0] + "\n")
    cov[vals[0]].append(float(vals[3]))
genes={}
for line in bedFile:
    vals=line.split()
    chrom=vals[0]
    start= int(int(vals[1])/100);
    end=int(int(vals[2])/100);
    gene=vals[3]
    if gene not in genes:
        genes[gene]=[0,0]
    
    if chrom not in cov:
        sys.stdout.write("0\n")
    else:
        total=sum(cov[chrom][start:end])
        genes[gene][0]+=total
        genes[gene][1]+= end-start

for gene in genes:
    if (genes[gene][1] > 0):
        sys.stdout.write(gene + "\t{:2.2f}\n".format(genes[gene][0] / genes[gene][1]))
#        sys.stdout.write("{:2.2f}\n".format(total/(meanCov*(end-start))))


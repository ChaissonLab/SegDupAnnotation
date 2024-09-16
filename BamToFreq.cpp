#include "htslib/hts.h"
#include "htslib/sam.h"
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include "htslib/faidx.h"
using namespace std;
int GetBest(int na, int nc, int ng, int nt, int &first) {

	return 0;
}
vector<int> NucMap;
int GetRefAndAlt(char refNuc, vector<int> &counts, int &ref, int &alt) {
  vector<int> sCounts=counts;
  sort(sCounts.begin(), sCounts.end());
  
  int second = sCounts[2];
  int first  = sCounts[3];
  int firstIndex=-1, secondIndex=-1;
  for (int i=0; i < 4; i++) {
    if (firstIndex == -1 and counts[i] == first) {
      firstIndex=i;
    }
    else if (secondIndex == -1 and counts[i] == second) {
      secondIndex=i;
    }
  }
  if (firstIndex == -1 or secondIndex == -1) {
    ref=0; alt=0;
    return 4;
  }
  int refNucIndex=NucMap[refNuc];
  if (first == 0 or second/first < 0.01) {
    ref=0;
    alt=0;
    return 4;
  }
  if (firstIndex == refNucIndex) {
    ref=first;
    alt=second;
    return secondIndex;
  }
  else {
    ref=second;
    alt=first;
    return firstIndex;
  }  
}

const char* nucs = "ACGTN";


int main(int argc, const char* argv[]) {
	if (argc < 3) {
		cout << "usage: bamToFreq file.bam regions.txt reference.fa" << endl
				 << "    The per nucleotide frequency will be calculated " << endl
				 << "    for the bam file for every region specified in regions.txt" << endl
				 << "    This file has one region per line, in the format chrom:start-end" << endl;
		exit(1);
	}
	NucMap.resize(256,4);
	NucMap[(int)'A']=0;
	NucMap[(int)'C']=1;
	NucMap[(int)'G']=2;
	NucMap[(int)'T']=3;	
	
	string bamFileName=argv[1];
	string regionFileName=argv[2];
	string referenceName=argv[3];
	//
	// Get the header of the bam file
	//
	htsFile *htsfp;

	htsfp = hts_open(bamFileName.c_str(),"r");
	bam_hdr_t *samHeader;			
	samHeader = sam_hdr_read(htsfp);

	faidx_t *fai = fai_load_format(referenceName.c_str(), FAI_FASTA);

	// 
	// Read index for random io
	//
	hts_idx_t *idx;
	if ((idx = sam_index_load(htsfp, bamFileName.c_str())) == 0) {
		cerr << "ERROR reading index" << endl;
		exit(0);
	}

	const htsFormat *fmt = hts_get_format(htsfp);
	if (fmt == NULL or (fmt->format != sam and fmt->format != bam)) {
		cout << "Cannot determine format of input reads." << endl;
		exit(1);
	}



	ifstream regionFile(regionFileName.c_str());

	string region;

  while ( (regionFile >> region) ) {
		cerr << "Parsing " << region << endl;
		hts_itr_t *regionIter = sam_itr_querys(idx, samHeader, region.c_str());
		int start;
		int end;
		string chrom;
		const char* chromEnd = hts_parse_reg(region.c_str(), &start, &end);
		start+=1;
		chrom=region.substr(0,chromEnd-region.c_str());
		bam1_t *b = bam_init1();
		int nReads=0;
		int regionLength= end-start;
		if (regionLength == 0) {
			continue;
		}

		vector<int> nA(regionLength, 0), nC(regionLength, 0), nT(regionLength, 0), nG(regionLength,0), nDel(regionLength, 0);
		while (bam_itr_next(htsfp, regionIter, b) > 0) {
			nReads+=1;
			if (nReads % 10000 == 0) { 
				cerr << "... " << nReads << endl;
			}
			int readLength = b->core.l_qseq;			
			if (readLength < 100) {
				continue;
			}
			if (b->core.qual < 10) {
				continue;
			}
			if (b->core.flag & 0x800) {
			  continue;
			}
			// Extract sequence

			vector<char> seq(readLength);
			uint8_t *q = bam_get_seq(b);
			for (int i=0; i < readLength; i++) {seq[i]=seq_nt16_str[bam_seqi(q,i)];	}
			uint32_t* cigar = bam_get_cigar(b);
			int refLen = bam_cigar2rlen(b->core.n_cigar, cigar);
			//			cout << bam_get_qname(b)  << "\t" << refLen << "\t" << (int) b->core.qual << "\t" << (int) b->core.flag << endl;			
			int qPos=0;
			int refPos = b->core.pos;
			int ci;
			int regionOffset=0;
			bool first=true;
			if (refPos >= start ) {
				regionOffset = refPos-start;
			}
			for (ci=0; ci < b->core.n_cigar; ci++) {
				int opLen=bam_cigar_oplen(cigar[ci]);
				int op=bam_cigar_op(cigar[ci]);

				if (op == BAM_CSOFT_CLIP) {
					qPos += opLen;
					continue;
				}
				if (op == BAM_CINS) {
					qPos += opLen;
					continue;
				}
				if (op == BAM_CDEL) {
					int stop=refPos+opLen;
					for (; refPos < end and refPos < stop; refPos++) {
						if (refPos >= start) {
							nDel[regionOffset]+=1;
							regionOffset++;
						}
					}
					continue;
				}
				if (op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF) {
					if (refPos + opLen <= start) {
						refPos += opLen;
						qPos += opLen;
						continue;
					}
					else {
						for (int p=0; p < opLen; p++) {
							if (refPos >= end){ 
								break;
							}
							if (refPos >= start) {
								first=false;
								char nuc=toupper(seq[qPos]);

								assert(regionOffset < nA.size());
								if (nuc == 'A') { nA[regionOffset]++;}
								if (nuc == 'C') { nC[regionOffset]++;}
								if (nuc == 'G') { nG[regionOffset]++;}
								if (nuc == 'T') { nT[regionOffset]++;}
								regionOffset++;
							}
							refPos++;
							qPos++;
						}
					}
				}
			}
			bam_destroy1(b);
			b = bam_init1();						
		}
		int chromSeqLen=0;
		char *chromSeq=faidx_fetch_seq(fai, chrom.c_str(), start, end, &chromSeqLen);
		vector<int> counts(4);
		for (int i = 0; i < nA.size(); i++) {
			int  len;
			//			char *nuc = faidx_fetch_seq(fai, chrom.c_str(), start+i,start+i, &len);
			int svn, sv;
			char refNuc=toupper(chromSeq[i]);
			int refCount, altCount;
			counts[0] = nA[i]; counts[1]=nC[i]; counts[2]=nG[i]; counts[3]=nT[i];
			svn= GetRefAndAlt(refNuc, counts, refCount, altCount);
			cout << chrom << "\t" << start+i << "\t" << nA[i] << "\t" << nC[i] << "\t" << nG[i] << "\t" << nT[i] << "\t" << nDel[i] << "\t0\t" << (char) toupper(chromSeq[i]) << "\t" << nucs[svn] << "\t" << refCount << "\t" << altCount << endl;
		}
		free(chromSeq);
	}
}

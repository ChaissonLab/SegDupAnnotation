#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "htslib/sam.h"
#include <cmath>
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << "Usage: gemBamToMask input.bam" << endl;
		exit(1);
	}
	string filename = argv[1];

	htsFile *htsfp;
	bam_hdr_t *samHeader;			


	htsfp = hts_open(filename.c_str(),"r");
	const htsFormat *fmt = hts_get_format(htsfp);


	samHeader = sam_hdr_read(htsfp);

	cerr << "starting " << endl;
	
	bam1_t *b = bam_init1();
	int res=1;
	res= sam_read1(htsfp, samHeader, b);
	long readIndex=0;
	long nCounted=0;
	cerr << "res: "<< res<< endl;
	while (res > 0) {
		readIndex+=1;
		long alnPos = b->core.pos;
		int tid=b->core.tid;
		if (alnPos >= 0) {
			nCounted+=1;
			uint8_t *xaData = bam_aux_get(b, "XA");
			if (xaData != 0) {
				char *xaString=bam_aux2Z(xaData);
				stringstream auxStrm((char*)xaString);
				string aln;
				int index=0;
				int nAux=0;
				while(std::getline(auxStrm, aln, ';')) { 				
					nAux+=1;
					if (nAux > 10000) {
						cerr <<"Woah nelly: " << readIndex << endl;
						exit(0);
					}
				}
				if (nAux>= 20) {
					cout << samHeader->target_name[tid] << "\t" << alnPos << "\t" << alnPos + b->core.l_qseq << "\t" << nAux << endl;
				}
			}
		}
		res = sam_read1(htsfp, samHeader, b); 
		if (readIndex %100000 == 0) {
			cerr << "proc " << readIndex /100000 << "M\t" << nCounted << endl;
		}
	}
}

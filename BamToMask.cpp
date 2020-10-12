#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "htslib/sam.h"
#include <cmath>

using namespace std;

int main(int argc, char* argv[]) {
	string filename = argv[1];

	htsFile *htsfp;
	bam_hdr_t *samHeader;			


	htsfp = hts_open(filename.c_str(),"r");
	const htsFormat *fmt = hts_get_format(htsfp);


	samHeader = sam_hdr_read(htsfp);

	
	vector<vector<int> > covBins;
	for (int i =0; i < samHeader->n_targets; i++) {
		covBins.push_back(vector<int>() );
		int last=covBins.size()-1;
		covBins[last].resize(samHeader->target_len[i]/100+1);
	}

	
	bam1_t *b = bam_init1();
	int res=1;
	res= sam_read1(htsfp, samHeader, b);
	long readIndex=0;
	long nCounted=0;
	while (res > 0) {
		readIndex+=1;
		long alnPos = b->core.pos;
		int tid=b->core.tid;
		if (alnPos >= 0) {
			covBins[tid][alnPos/100]+=1;
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
				}
				if (nAux>= 20) {
					cout << samHeader->target_name[tid] << "\t" << alnPos << "\t" << alnPos + b->core.l_qseq << "\t" << nAux << endl;
				}
			}
			res = sam_read1(htsfp, samHeader, b); 
			if (readIndex %1000000 == 0) {
				cerr << "proc " << readIndex /1000000 << "M\t" << nCounted << endl;
			}
		}
	}
}

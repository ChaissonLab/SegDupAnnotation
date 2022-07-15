#include <cinttypes>
#include <assert.h>
#include "htslib/htslib/kseq.h"
#include <sstream>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread);

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

using namespace std;
int main(int argc, char* argv[]) {
	if (argc == 1) {
		cout << "usage: mask2bed input.fasta output.bed" << endl;
		exit(1);
	}
	string outFileName=argv[2];
	ifstream testFile(outFileName.c_str());
	if (testFile.good()) {
		cout << "ERROR, output file " << outFileName << " already exists." << endl;
		exit(1);
	}
	ofstream outFile(argv[2]);

	gzFile fastaFile;
	kseq_t * ks;
	map<string, char*> nameToSeq;
	map<string, int>   nameToLen;
	cerr << argv[1] << endl;
	fastaFile = gzopen(argv[1], "r");
	ks = kseq_init(fastaFile);
	vector<string> nameOrder;
	int pi;
	while (kseq_read(ks) >= 0) {
		cerr << ks->name.s<< endl;
		int i=0;
		int j=0;
		while (i < ks->seq.l) {
			pi=i;
			while (i < ks->seq.l and ks->seq.s[i] >= 'A' and ks->seq.s[i] <= 'Z') {
				i+=1;
			}
			j=i;
			while (j < ks->seq.l and ks->seq.s[j] >= 'a' and ks->seq.s[j] <= 'z') {
				j+=1;
			}
			if (j > i) {
				outFile << ks->name.s << "\t" << i << "\t" << j << endl;
			}
			i=j;
			// Skip past strange characters.
			if (pi == i) {
				i++;
			}
		}
	}
}

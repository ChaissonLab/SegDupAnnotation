#include <cinttypes>

#include "htslib/kseq.h"

#include <zlib.h>

KSEQ_INIT(gzFile, gzread);

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
int main(int argc, char* argv[]) {
	if (argc == 1) {
		cout << "usage: comask output.fasta input1.fa input2.fa ..." << endl;
		exit(1);
	}

	string outFileName = argv[1];
	ifstream testFile(outFileName.c_str());
	if (testFile.good()) {
		cout << "ERROR, output file " << outFileName << " already exists." << endl;
		exit(1);
	}
	ofstream outFile(outFileName.c_str());
	vector< string > fileNames;
	vector<	gzFile > fastaFiles;
	vector<	kseq_t *> ks;

	for (int i = 2; i < argc; i++) {
		fileNames.push_back(argv[i]);
		fastaFiles.push_back(gzopen(argv[i], "r"));
		ks.push_back(kseq_init(fastaFiles[fastaFiles.size()-1]));
	}

	if (ks.size() == 0) {
		cout << "No input files " << endl;
		exit(1);
	}
	while (kseq_read(ks[0]) >= 0) {
		cerr << ks[0]->name.s<< endl;
		for (int j = 1; j < ks.size(); j++) {
			if (kseq_read(ks[j]) == 0) {
				cout << "ERROR, could not fetch sequence " << ks[0]->name.s << " from " << fileNames[j] << endl;
				exit(0);
			}
			if (strcmp(ks[j]->name.s, ks[0]->name.s) != 0) {
				cout << "ERROR, input file " << fileNames[j] << " is out of sync with " << fileNames[0] << endl;
				exit(0);
			}
		}
		for (int p=0; p < ks[0]->seq.l; p ++) {
			for (int j=1; j < ks.size(); j++) {
				if (ks[j]->seq.s[p] >= 'a' && ks[j]->seq.s[p] <= 'z') {
					ks[0]->seq.s[p] = ks[j]->seq.s[p];
				}
			}
		}
		outFile << ">" << ks[0]->name.s << endl;
		long p=0;
		while (p < ks[0]->seq.l) {
			long int end =  min(p+60, (long int) ks[0]->seq.l);
			string sub(&ks[0]->seq.s[p], end-p);
			outFile << sub << endl;
			p=end;
		}
	}

}

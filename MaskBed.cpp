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
		cout << "usage: bemask input.fasta mask.bed output.fasta" << endl;
		exit(1);
	}

	string outFileName = argv[3];
	ifstream testFile(outFileName.c_str());
	if (testFile.good()) {
		cout << "ERROR, output file " << outFileName << " already exists." << endl;
		exit(1);
	}
	ofstream outFile(outFileName.c_str());
	ifstream bedFile(argv[2]);
	gzFile fastaFile;
	kseq_t * ks;
	map<string, char*> nameToSeq;
	map<string, int>   nameToLen;
	fastaFile = gzopen(argv[1], "r");
	ks = kseq_init(fastaFile);
	vector<string> nameOrder;
	while (kseq_read(ks) >= 0) {
		cerr << ks->name.s<< endl;
		nameToSeq[ks->name.s] = new char[ks->seq.l];
		memcpy(nameToSeq[ks->name.s], ks->seq.s, ks->seq.l);
		nameToLen[ks->name.s] = ks->seq.l;
		nameOrder.push_back(ks->name.s);
	}
	int lineNumber=0;
	while (bedFile) {
		
		string line;
		getline(bedFile, line);
		stringstream strm(line);
		string chrom;
		int start, end;
		strm >> chrom >> start >> end;
		if (chrom == "") {
			break;
		}
		if (nameToSeq.find(chrom) == nameToSeq.end()) {
			cerr << "ERROR! " << chrom << " not found " << endl;
			exit(1);
		}
		char *seq=nameToSeq[chrom];
		assert(end <= nameToLen[chrom]);
		int  l=nameToLen[chrom];
		for (int i=start; i < end; i++) {
			seq[i] = tolower(seq[i]);
		}
		++lineNumber;
	}

	for (int s=0; s < nameOrder.size(); s++) {
		char *seq=nameToSeq[nameOrder[s]];
		int  l=nameToLen[nameOrder[s]];
		outFile << ">" << nameOrder[s] << endl;
		int p=0;
		while (p < l) {
			int end = min(p+60, (int) l);
			string sub(&seq[p], end-p);
			outFile << sub << endl;
			p=end;
		}
	}

}

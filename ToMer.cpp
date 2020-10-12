#include <iostream>
#include <string>
#include <vector>

using namespace std;
int main( int argc, char* argv[]) {
	if (argc < 2) {
		cout << "usage: tomer k  < infile.fastq" << endl;
		exit(0);
	}
	int kmer=atoi(argv[1]);

	while (cin) {
		string title, seq,sep,qual;
		getline(cin,title);
		getline(cin,seq);
		getline(cin,sep);
		getline(cin,qual);
		if (seq.size() >= kmer) {
			for (int i=0; i < seq.size()-kmer + 1; i++) {
				string ks=seq.substr(i,kmer);
				cout << ">" << title << "/" << i << endl;
				cout << ks << endl;
			}
		}
	}
}

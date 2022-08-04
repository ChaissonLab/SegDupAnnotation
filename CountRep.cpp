#include <iostream>
#include <string>
#include <vector>
using namespace std;

// purpose: Calculate percent assembly masked.
// method: Report number of lower and upper case chars in input.
// input: (assembly) fasta file
// output: <fasta entry header> <(num lowercase chars)/(num total chars)> <num lowercase chars in entry> <num uppercase chars in entry>

int main() {
	long int nl=0, nh=0, other=0;

	string prevTitle = "";
	while (cin) {
		string line;
		cin >> line;
		if (line.size() > 0) {
			if (line[0] == '>') {
				if (prevTitle != "") {
					int s=nl+nh;
					if (s > 0) {
						cout << prevTitle << "\t"<< nl / ((float)s) << "\t" << nl << "\t" << nh << endl;
					}
					else {
						cout << prevTitle << "\t0\t0\t0" << endl;
					}
				}
				nl=0; nh=0;
				prevTitle = line;
			} else {
				for (int i=0; i < line.size(); i++) {
					if (line[i] >= 'a' && line[i] <= 'z') {
						nl+=1;
					}
					else {
						nh+=1;
					}
				}
			}
		}
	}
	int s=nl+nh;
	if (s > 0) {
		cout << prevTitle << "\t"<< nl / ((float)s) << "\t" << nl << "\t" << nh<< endl;
	}
	else {
		cout << prevTitle << "\t0\t0\t0" << endl;
	}
}
		

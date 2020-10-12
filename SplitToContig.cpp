#include <iostream>
#include <string>
#include <vector>
using namespace std;
int main() {
	int nl=0, nh=0, other=0;
	vector<int> low(256,0), high(256,0);
	string prevTitle = "";
	int start, cur;
	char prev='Z';
	string contig;
	while (cin) {
		string line;
		cin >> line;
		if (line.size() > 0 and line[0] == '>') {
			cerr << "proc " << line << endl;
			if (prev == 'Z') {
				cout << cur-prev << endl;
			}
			cur=0;
			start=0;
			prev='A';
		}
		else {
			for (int i=0; i < line.size(); i++) {
				if ( line[i] == 'N' or line[i] == 'n') {
					if (prev != 'n') {
						cout << cur - start<< endl;
					}
					prev='n';
				}
				else {
					cur++;
					if (prev == 'n') {						
						start = 0;
						cur=0;
					}
					prev=line[i];
				}
			}
		}
	}
}
		

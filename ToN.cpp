#include <iostream>
#include <string>
#include <vector>
using namespace std;
int main() {
	int nl=0, nh=0, other=0;
	vector<int> low(256,0), high(256,0);
	string prevTitle = "";
	while (cin) {
		string line;
		cin >> line;
		if (line.size() > 0 and line[0] == '>') {
			cout << line << endl;
		}
		else {
			for (int i=0; i < line.size(); i++) {
				if (line[i] != 'A' or line[i] != 'C' or line[i] != 'G' or line[i] != 'T') {
					line[i] =  'N';
				}
			}
			cout << line << endl;
		}
	}
}

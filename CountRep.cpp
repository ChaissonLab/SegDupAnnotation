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
		}
		for (int i=0; i < line.size(); i++) {
			if (line[i] >= 'a' && line[i] <= 'z') {
				nl+=1;
			}
			else {
				nh+=1;
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
		

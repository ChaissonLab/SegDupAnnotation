#include <cinttypes>
#include <assert.h>

#include <sstream>




#include <vector>
#include <string>
#include <iostream>
#include <fstream>


using namespace std;
int main(int argc, char* argv[]) {
  
  if (argc == 1) {
    cout << "usage: upper input.fasta output.fasta" << endl;
    exit(1);
  }
  string inFileName  = argv[1];
  string outFileName = argv[2];
  ifstream testFile(outFileName.c_str());
  if (testFile.good()) {
    cout << "ERROR, output file " << outFileName << " already exists." << endl;
    exit(1);
  }
  ofstream outFile(outFileName.c_str());

  ifstream fastaFile(inFileName.c_str());

   
  vector<string> nameOrder;
  string line;

  while (getline(fastaFile, line)) {
    if (line.size() > 0 and line[0] == '>') {
      outFile << line << endl; 
    }
    else {
      for (int i=0; i < line.size(); i++) {
	char c=line[i];
	if (c == 'a' || c == 'c' || c == 'g' || c=='t') {
	  line[i] = 'N';
	}
      }
      outFile << line << endl;
    }
  }
}

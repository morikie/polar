#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>
#include "fastaReader.hpp"
#include "utr3MutationFinder.hpp"
#include "polar.hpp"



using namespace std;

int main (int args, char * argv[]) {
	
	fs::path f = "UTR3_sequences_test.txt";
        FastaReader newReader(f);

//	ifstream myFile("UTR3_sequences_test.txt");
//	string line;
//	vector<string> myLines;
//	while (getline(myFile, line)) {
//		myLines.push_back(line);
//	}
//	for (size_t i = 0; i < 2; ++i) {
//		cerr << myLines[i] << endl;
//	}
	return EXIT_SUCCESS;
}

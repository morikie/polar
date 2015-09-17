#include <iostream>
#include <boost/filesystem.hpp>
#include "../src/utr3Finder.hpp"
#include "readKnownPolyA.hpp"
#include "acc_test.hpp"

namespace fs = boost::filesystem;


int main (int argc, char * argv[]) {
	fs::path knownPolyA = "../perf_testing/knownPolyAtranscript.txt";
	std::vector<KnownPolyA> kPolyAvec;
	if (fs::exists(knownPolyA)) {
		readKnownPolyA(knownPolyA, kPolyAvec);
	} else {
		std::cerr << "error: File not found." << std::endl;
		return EXIT_FAILURE;
	}
	std::cerr << kPolyAvec.size() << std::endl;
	/* for (std::vector<KnownPolyA>::iterator it = kPolyAvec.begin(); it != kPolyAvec.end(); it++) {
		
		std::cerr << it->id << ", Poly(A) position(s): ";
		for (size_t i = 0; i < it->polyApos.size(); i++) {
			std::cerr << it->polyApos[i] << ", ";
		}
		std::cerr << std::endl;
	}
	*/
	std::vector<SeqStruct> seqStt;

	buildSeqStruct(seqStt, kPolyAvec);

}

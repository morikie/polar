#include <iostream>
#include <boost/filesystem.hpp>
#include "../src/hgvsParser.hpp"
#include "../src/utr3Finder.hpp"
#include "../src/utr3FinderFuzzy.hpp"
#include "../src/utr3FinderNaive.hpp"
#include "readKnownPolyA.hpp"
#include "acc_test.hpp"

namespace fs = boost::filesystem;


int main (int argc, char * argv[]) {
	fs::path knownPolyA = "../perf_testing/knownPolyAtranscript.txt";
	std::vector<SeqStruct> seqStt;
	std::vector<KnownPolyA> kPolyAvec;
	std::vector<Utr3Finder*> u3F;
	if (fs::exists(knownPolyA)) {
		readKnownPolyA(knownPolyA, kPolyAvec);
	} else {
		std::cerr << "error: File not found." << std::endl;
		return EXIT_FAILURE;
	}
	/* for (std::vector<KnownPolyA>::iterator it = kPolyAvec.begin(); it != kPolyAvec.end(); it++) {
		
		std::cerr << it->id << ", Poly(A) position(s): ";
		for (size_t i = 0; i < it->polyApos.size(); i++) {
			std::cerr << it->polyApos[i] << ", ";
		}
		std::cerr << std::endl;
	}
	*/
	buildSeqStruct(seqStt, kPolyAvec);
	
	for (auto it = seqStt.begin(); it != seqStt.end(); ++it) {
		u3F.push_back(new Utr3FinderNaive(*it));
	}
	size_t correctPredict = 0;
	size_t total = u3F.size();
	for (size_t i = 0; i < u3F.size(); i++) {
		size_t predictedPos = u3F[i]->getPolyaMotifPos()[0].pos;
		for (size_t j = 0; j < kPolyAvec[i].polyApos.size(); j++) {
			if (predictedPos == kPolyAvec[i].polyApos[j]) {
				correctPredict++;
				break;
			}
			if (j + 1 == kPolyAvec[i].polyApos.size()) {
				std::cerr << "reference polyA: " << kPolyAvec[i].polyApos[j] << ", ";
				u3F[i]->writeInfo();
			}
		}
	}
	
	std::cerr << "total: "<< total << ", correct predictions: " << correctPredict << std::endl;
	double corrPerc = static_cast<double>(correctPredict) / static_cast<double>(total);
	std::cerr << "relative correctness: " << corrPerc << std::endl;
}


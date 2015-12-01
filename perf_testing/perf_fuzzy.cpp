#include <iostream>
#include <boost/filesystem.hpp>
#include "../src/hgvsParser.hpp"
#include "../src/refGeneParser.hpp"
#include "../src/utr3Finder.hpp"
#include "../src/utr3FinderFuzzy.hpp"
#include "../src/utr3FinderNaive.hpp"
#include "readKnownPolyA.hpp"
#include "perf_fuzzy.hpp"

namespace fs = boost::filesystem;


int main (int argc, char * argv[]) {
	fs::path refGeneFile = "ucsc_data/refGene.txt";
	RefGeneParser test(refGeneFile);

	RefGeneProperties testProp = test.getValueByKey("NM_001012993");

	if (testProp.chr != "chr9") std::cerr << "bad chromosome" << std::endl;
	if (testProp.txStart != 112961845) std::cerr << "bad tx start" << std::endl;
	if (testProp.name2 != "C9orf152") std::cerr << "bad name2" << std::endl;
	if (testProp.cdsEndStat != "cmpl") std::cerr << "bad status" << std::endl;

	std::cerr << "test done." << std::endl;

}


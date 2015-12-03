#include <iostream>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/spirit/include/qi.hpp>
#include "../src/hgvsParser.hpp"
#include "../src/refGeneParser.hpp"
#include "../src/utr3Finder.hpp"
#include "../src/utr3FinderFuzzy.hpp"
#include "../src/utr3FinderNaive.hpp"
#include "readKnownPolyA.hpp"
#include "perf_fuzzy.hpp"

namespace fs = boost::filesystem;
namespace qi = boost::spirit::qi;


int main (int argc, char * argv[]) {
	fs::path refGeneFile = "ucsc_data/refGene.txt";
	fs::path knownPolyA = "../perf_testing/knownPolyAtranscript.txt";
	std::vector<KnownPolyA> knownPolyAvec;
	std::map<std::string, std::vector<size_t> > truePositives;
	RefGeneParser refGene(refGeneFile);
	
	readKnownPolyA(knownPolyA, knownPolyAvec);

	BOOST_FOREACH (KnownPolyA & knownPolyA, knownPolyAvec) {
		std::string baseSeqId;
		qi::parse(knownPolyA.id.begin(), knownPolyA.id.end(), *~qi::char_('.'), baseSeqId);
		size_t txLength = knownPolyA.seq.size();
		RefGeneProperties refGeneProp = refGene.getValueByKey(baseSeqId);
		if (knownPolyA.seq == std::string()) continue;
		BOOST_FOREACH (size_t & pos, knownPolyA.polyApos) {
			size_t geneticPos;
			if (refGeneProp.strand == "+") {
				geneticPos = refGeneProp.txEnd - (txLength - pos);
			} else {
				geneticPos = refGeneProp.txStart - (txLength - pos);
			}
			truePositives[refGeneProp.chr].push_back(geneticPos);
		}
	}

	for (auto iter = truePositives.begin(); iter != truePositives.end(); iter++) {
		std::cerr << iter->first << ": ";
		for (auto it = iter->second.begin(); it != iter->second.end(); it++) {
			std::cerr << *it << ", ";
		}
		std::cerr << std::endl;
	}
}


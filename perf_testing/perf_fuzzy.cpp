#include <climits>
#include <iostream>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/spirit/include/qi.hpp>
#include <seqan/seq_io.h>
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
	typedef std::string key;
	typedef std::string strand;
	typedef size_t position;
	typedef std::pair<position, strand> posStrandPair;
	std::map<key, std::vector<posStrandPair> > truePositives;
	fs::path refGeneFile = "ucsc_data/refGene.txt";
	fs::path knownPolyA = "../perf_testing/knownPolyAtranscript.txt";
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	std::vector<KnownPolyA> knownPolyAvec;
	RefGeneParser refGene(refGeneFile);	
	size_t correctMatches = 0;
	size_t incorrectMatches = 0;
	size_t unknMotives = 0;
	size_t analyzedSequences = 0;

	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		std::cerr << "could not open index file for " << referenceGenome << std::endl;
	}
	readKnownPolyA(knownPolyA, knownPolyAvec);

	BOOST_FOREACH (KnownPolyA & knownPolyA, knownPolyAvec) {
		std::string baseSeqId;
		qi::parse(knownPolyA.id.begin(), knownPolyA.id.end(), *~qi::char_('.'), baseSeqId);
		RefGeneProperties refGeneProp = refGene.getValueByKey(baseSeqId);
		if (knownPolyA.seq.empty() || refGeneProp.chr.empty()) continue;
		if (refGeneProp.chr.size() > 6) continue;
		
		size_t txLength = 0;
		for (size_t i = 0; i < refGeneProp.exonStarts.size(); i++) {
			txLength += refGeneProp.exonEnds[i] - refGeneProp.exonStarts[i];
		}

		BOOST_FOREACH (size_t & pos, knownPolyA.polyApos) {
			size_t geneticPos;
			if (refGeneProp.strand == "+") {
				geneticPos = refGeneProp.txEnd - (txLength - pos);
			} else {
				geneticPos = refGeneProp.txStart + (txLength - pos - 6);
			}
			truePositives[refGeneProp.chr].push_back(std::make_pair(geneticPos, refGeneProp.strand));
		}
	}
	
	for (auto it = truePositives.begin(); it != truePositives.end(); it++) {
		for (auto vecIt = it->second.begin(); vecIt != it->second.end(); vecIt++) {
			size_t & pos = vecIt->first;
			std::string & strand = vecIt->second;
			seqan::CharString temp;
			unsigned idx = UINT_MAX;
			if (it->first == "chrX") {
				idx = 22;
			} else if (it->first == "chrY") {
				idx = 23;
			} else {
				if (! qi::parse(it->first.begin(), it->first.end(), qi::omit[*qi::alpha] >> qi::uint_, idx)) {
					throw std::invalid_argument("error parsing chromosome value from vcf");
				}
				//adjusting idx to 0-starting map
				idx--;
			}
			
			seqan::readRegion(temp, faiIndex, idx, pos - 100, pos + 100);
			SeqStruct ss = {
				std::string(seqan::toCString(temp)),
				boost::none,
				boost::none,
				boost::none,
				boost::none,
				boost::none,
				boost::none,
				boost::none
			};
			//std::cerr << "seq.size(): " << ss.seq.size() << "; " << ss.seq << std::endl;
			Utr3FinderFuzzy u3Fuzzy(ss);
			std::vector<Utr3Finder::Utr3FinderResult> u3FuzzyResVector = u3Fuzzy.getPolyaMotifPos();
			


			for (auto resultIt = u3FuzzyResVector.begin(); resultIt != u3FuzzyResVector.end(); resultIt++) {
				if (resultIt->pos == 100 && resultIt->strand == "+" && strand == "+") {
					correctMatches++;
				} else if (resultIt->pos == 94 && resultIt->strand == "-" && strand == "-") {
					correctMatches++;
				} else {
					incorrectMatches++;
					std::string motif = ss.seq.substr(100, 6);
					if (strand == "-") {
						std::string temp;
						for (auto it = motif.rbegin(); it != motif.rend(); it++) {
							temp.push_back(Utr3Finder::complement(*it));
						}
						motif = temp;
					}
					
					std::cerr << resultIt->truthValue << "  " << motif << std::endl;
					if (Utr3FinderFuzzy::dseLocMap.find(motif) == Utr3FinderFuzzy::dseLocMap.end()) {
						unknMotives++;			
					}
				}
			}
			analyzedSequences++;
		}
	}
	std::cerr << "correct predictions: " << correctMatches << std::endl;
	std::cerr << "incorrect predictions: " << incorrectMatches << std::endl;
	std::cerr << "unknown motives: " << unknMotives << std::endl;
	std::cerr << "total sequences analyzed (true positives): " << analyzedSequences << std::endl; 
}


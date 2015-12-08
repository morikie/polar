#include <algorithm>
#include <climits>
#include <iostream>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/spirit/include/qi.hpp>
#include <seqan/seq_io.h>
#include "../src/hgvsParser.hpp"
#include "../src/polarUtility.hpp"
#include "../src/refGeneParser.hpp"
#include "../src/utr3Finder.hpp"
#include "../src/utr3FinderFuzzy.hpp"
#include "../src/utr3FinderNaive.hpp"
#include "readKnownPolyA.hpp"
#include "perf_fuzzy.hpp"

namespace fs = boost::filesystem;
namespace qi = boost::spirit::qi;

using namespace polar::utility;

int main (int argc, char * argv[]) {
	typedef std::string key;
	typedef std::string strand;
	typedef size_t position;
	typedef std::pair<position, strand> posStrandPair;
	//map that stores the position and strand for each chromosome (key)
	std::map<key, std::vector<posStrandPair> > truePositives;
	fs::path refGeneFile = "ucsc_data/refGene.txt";
	fs::path knownPolyA = "../perf_testing/knownPolyAtranscript.txt";
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	std::vector<KnownPolyA> knownPolyAvec;
	RefGeneParser refGene(refGeneFile);	
	
	size_t numTruePositives = 0; //found true positives/negatives
	size_t unknMotives = 0; //motives that aren't searched for or refGene.txt contains faulty mappings
	size_t totalTruePositives = 0; //total number of true positives

	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		std::cerr << "could not open index file for " << referenceGenome << std::endl;
	}
	readKnownPolyA(knownPolyA, knownPolyAvec);
	
	//filling up the true positives map
	BOOST_FOREACH (KnownPolyA & knownPolyA, knownPolyAvec) {
		std::string baseSeqId;
		qi::parse(knownPolyA.id.begin(), knownPolyA.id.end(), *~qi::char_('.'), baseSeqId);
		RefGeneProperties refGeneProp = refGene.getValueByKey(baseSeqId);
		//ignoring empty entries
		if (knownPolyA.seq.empty() || refGeneProp.chr.empty()) continue;
		//ignoring entries from not processable chromosome descriptions (e.g. chr6_hash_map... or similar)
		if (refGeneProp.chr.size() > 6) continue;
		
		//calculate size of the transcript by using exon start/end values (using sequence length doesn't work)
		size_t txLength = 0;
		for (size_t i = 0; i < refGeneProp.exonStarts.size(); i++) {
			txLength += refGeneProp.exonEnds[i] - refGeneProp.exonStarts[i];
		}
		
		//mapping transcript position of the PAS to the genomic position (each position is a true positive)
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
	//evaluating every true positive
	for (auto it = truePositives.begin(); it != truePositives.end(); it++) {
		//iterating over every true PAS
		for (auto vecIt = it->second.begin(); vecIt != it->second.end(); vecIt++) {
			size_t & pos = vecIt->first;
			std::string & strand = vecIt->second;
			seqan::CharString temp;
			size_t idx = polar::utility::getFastaIndex(it->first);
			//copy the genomic sequence 100 bases around the genomic position of the PAS (200nt long)
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
			//evaluating the sequence
			Utr3FinderFuzzy u3Fuzzy(ss);
			std::vector<Utr3Finder::Utr3FinderResult> u3FuzzyResVector = u3Fuzzy.getPolyaMotifPos();
			//analyzing the results
			for (auto resultIt = u3FuzzyResVector.begin(); resultIt != u3FuzzyResVector.end(); resultIt++) {
				if (resultIt->pos == 100 && resultIt->strand == "+" && strand == "+") {
					numTruePositives++;
				} else if (resultIt->pos == 94 && resultIt->strand == "-" && strand == "-") {
					numTruePositives++;
				} else {
					std::string motif = ss.seq.substr(100, 6);
					if (strand == "-") {
						std::string temp;
						for (auto it = motif.rbegin(); it != motif.rend(); it++) {
							temp.push_back(Utr3Finder::complement(*it));
						}
						motif = temp;
					}
					
					//std::cerr << resultIt->truthValue << "  " << motif << std::endl;
					if (Utr3FinderFuzzy::dseLocMap.find(motif) == Utr3FinderFuzzy::dseLocMap.end()) {
						unknMotives++;			
					}
				}
			}
			totalTruePositives++;
		}
	}
	std::cerr << "-------Sensitivity test-------" << std::endl;
	std::cerr << "correct predictions (found true positives): " << numTruePositives << std::endl;
	std::cerr << "incorrect predictions (false negatives): " << totalTruePositives - unknMotives - numTruePositives << std::endl;
	std::cerr << "total sequences analyzed (total true positives): " << totalTruePositives << std::endl; 
	std::cerr << "unknown motives: " << unknMotives << std::endl;
	double sensitivity = static_cast<double>(numTruePositives) / (totalTruePositives - unknMotives);
	std::cerr << "sensitivity (w/o unknown motives): " << sensitivity << std::endl;
	
	size_t numTrueNegatives = 0; //true negatives
	size_t numFalsePositives = 0; //"matches" found around a true positives/negatives  (not further analyzed atm)
	size_t totalTrueNegatives = 0; //total number of true positives

	std::vector<std::string> transcriptVector = refGene.getKeys();
	std::map<key, std::vector<posStrandPair> > trueNegatives;
	
	size_t hitCounter = 0;
	const size_t maxMatches = 30000;
	for (auto txIter = transcriptVector.begin(); txIter != transcriptVector.end(); txIter++) {
		RefGeneProperties refGeneProp = refGene.getValueByKey(*txIter);
		if (refGeneProp.chr.size() > 6 || refGeneProp.chr.empty()) continue;

		const std::vector<unsigned int> & txStartVector = refGeneProp.exonStarts;
		const std::vector<unsigned int> & txEndVector = refGeneProp.exonEnds;
		for (unsigned int i = 0; i < txStartVector.size(); i++) {
			seqan::CharString temp;
			size_t idx = polar::utility::getFastaIndex(refGeneProp.chr);	
			seqan::readRegion(temp, faiIndex, idx, txStartVector[i], txEndVector[i]);
			if (refGeneProp.strand == "-") seqan::reverseComplement(temp);
			std::string exonSequence(seqan::toCString(temp));

			auto hexamersIter = Utr3Finder::hexamers.begin();
			for (; hexamersIter != Utr3Finder::hexamers.end(); hexamersIter++) {
				auto match = exonSequence.begin();
				while ((match = std::search(match, exonSequence.end(),
					hexamersIter->begin(), hexamersIter->end())) != exonSequence.end()) {
					size_t geneticPos;
					if (refGeneProp.strand == "+") {
						geneticPos = std::distance(exonSequence.begin(), match) + txStartVector[i];
					} else {
						geneticPos = txEndVector[i] - std::distance(exonSequence.begin(), match);
					}
					trueNegatives[refGeneProp.chr].push_back(std::make_pair(geneticPos, refGeneProp.strand));
					match++;
					hitCounter++;
					if (hitCounter >= maxMatches) goto stopLoop;
				}
			}
		}
	}
	stopLoop:	
	//evaluating every true negative
	for (auto it = trueNegatives.begin(); it != trueNegatives.end(); it++) {
		//iterating over every true PAS
		for (auto vecIt = it->second.begin(); vecIt != it->second.end(); vecIt++) {
			size_t & pos = vecIt->first;
			std::string & strand = vecIt->second;
			seqan::CharString temp;
			//mapping chromosome to id number for use in the fasta index
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
			//copy the genomic sequence 100 bases around the genomic position of the putative PAS (200nt long)
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
			//evaluating the sequence
			Utr3FinderFuzzy u3Fuzzy(ss);
			std::vector<Utr3Finder::Utr3FinderResult> u3FuzzyResVector = u3Fuzzy.getPolyaMotifPos();

			std::cerr << ss.seq << std::endl;
			bool foundMatch = false;
			//analyzing the results
			for (auto resultIt = u3FuzzyResVector.begin(); resultIt != u3FuzzyResVector.end(); resultIt++) {
				if (resultIt->pos == 100 && resultIt->strand == "+" && strand == "+") {
					numFalsePositives++;
					foundMatch = true;
				} else if (resultIt->pos == 94 && resultIt->strand == "-" && strand == "-") {
					numFalsePositives++;
					foundMatch = true;	
				}
			}
			if (! foundMatch) numTrueNegatives++;
			totalTrueNegatives++;
		}
	}

	std::cerr << "-------Specifity test-------" << std::endl;
	std::cerr << "correct predictions (found true negatives): " << numTrueNegatives << std::endl;
	std::cerr << "incorrect predictions (false positives): " << numFalsePositives << std::endl;
	std::cerr << "total sequences analyzed (total true positives): " << totalTrueNegatives << std::endl; 
	//std::cerr << "incorrect predictions: " << numFalsePositives << std::endl;
	double specifity = static_cast<double>(numTrueNegatives) / totalTrueNegatives;
	std::cerr << "specifity: " << specifity << std::endl;
}


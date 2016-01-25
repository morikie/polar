#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/spirit/include/qi.hpp>
#include <seqan/seq_io.h>
#include "../src/hgvsParser.hpp"
#include "../src/polarUtility.hpp"
#include "../src/refGeneParser.hpp"
#include "../src/utr3Finder.hpp"
#include "../src/utr3FinderFuzzy.hpp"
#include "../src/utr3FinderNaive.hpp"
#include "createTNset.hpp"
#include "createTPset.hpp"
#include "readKnownPolyA.hpp"
#include "perf_fuzzy.hpp"

namespace fs = boost::filesystem;
namespace qi = boost::spirit::qi;


int main (int argc, char * argv[]) {
	typedef std::string key;
	typedef std::string strand;
	typedef size_t position;
	typedef std::pair<position, strand> posStrandPair;
	typedef std::pair<size_t, size_t> range;
	//map that stores position and strand of a match sorted by chromosome (key)
	std::map<key, std::vector<posStrandPair> > truePositives;
	fs::path refGeneFile = "ucsc_data/refGene.txt";
	fs::path knownPolyA = "../perf_testing/knownPolyAtranscript.txt";
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	fs::path ucscMappedTx = "ucsc_data/ucsc_txRefSeq.txt";
	fs::path tpFasta = "../perf_testing/tpSet.fa";
	fs::path tnFasta = "../perf_testing/tpSet.fa";
	std::unordered_map<std::string, size_t> ucscTxRefSeq = polar::utility::getTxRefSeqAccessions(ucscMappedTx);
	std::vector<KnownPolyA> knownPolyAvec;
	RefGeneParser refGene(refGeneFile);	
	std::cerr << "#accession numbers: " << ucscTxRefSeq.size() << std::endl;
	std::unordered_map<std::string, double> thresholdMap = {
		{std::string("aataaa"), 0.0},
		{std::string("attaaa"), 0.0},
		{std::string("tataaa"), 0.0},
		{std::string("agtaaa"), 0.0},
		{std::string("aagaaa"), 0.0},
		{std::string("aatata"), 0.0},
		{std::string("aataca"), 0.0},
		{std::string("cataaa"), 0.0},
		{std::string("gataaa"), 0.0},
		{std::string("aatgaa"), 0.0},
		{std::string("actaaa"), 0.0},
		{std::string("aataga"), 0.0}
	};
	
	if (! fs::exists(tpFasta)) {
		createTPset(tpFasta);
	}
	if (! fs::exists(tnFasta)) {
		createTNset(tnFasta);
	}

	std::vector<double> sensitivityVec;
	std::vector<double> specificityVec;
	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		std::cerr << "could not open index file for " << referenceGenome << std::endl;
	}
	readKnownPolyA(knownPolyA, knownPolyAvec);
	//creating true positive data set
	BOOST_FOREACH (KnownPolyA & knownPolyA, knownPolyAvec) {
		std::pair<std::string, size_t> txAndPatchPair;
		qi::parse(knownPolyA.id.begin(), knownPolyA.id.end(), *~qi::char_('.') >> '.' >> qi::uint_, txAndPatchPair);
		auto mapValue = ucscTxRefSeq.find(txAndPatchPair.first);
		if (mapValue == ucscTxRefSeq.end() || 
			mapValue->second != txAndPatchPair.second ||
			txAndPatchPair.first.find("NM_") == std::string::npos) {
			continue;
		}
		RefGeneProperties refGeneProp = refGene.getValueByKey(txAndPatchPair.first);
		size_t idx = polar::utility::getFastaIndex(refGeneProp.chr);
		//ignoring empty entries
		if (knownPolyA.seq.empty() || refGeneProp.chr.empty()) continue;
		//ignoring entries from not processable chromosome descriptions (e.g. chr6_hash_map etc.)
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
			
			if (refGeneProp.strand == "+") {
				seqan::CharString motifAtPos;
				seqan::readRegion(motifAtPos, faiIndex, idx, geneticPos, geneticPos + 6);
				std::string stdMotifAtPos(seqan::toCString(motifAtPos));
				auto findIt = thresholdMap.find(stdMotifAtPos);
				if (findIt == thresholdMap.end()) {
					continue;
				}
			} else {
				seqan::CharString motifAtPos;
				seqan::readRegion(motifAtPos, faiIndex, idx, geneticPos - 6, geneticPos);
				seqan::reverseComplement(motifAtPos);
				seqan::toLower(motifAtPos);
				std::string stdMotifAtPos(seqan::toCString(motifAtPos));
				auto findIt = thresholdMap.find(stdMotifAtPos);
				if (findIt == thresholdMap.end()) {
					continue;
				}
			}
			bool isDuplicate = false;
			auto mapVecIter = truePositives[refGeneProp.chr].begin();
			for(; mapVecIter != truePositives[refGeneProp.chr].end(); mapVecIter++) {
				if (mapVecIter->first == geneticPos) {
					isDuplicate = true;
					break;
				}
			}
			if (! isDuplicate) {
				truePositives[refGeneProp.chr].push_back(std::make_pair(geneticPos, refGeneProp.strand));
			}
		}
	}
	//creating true negative data set
	std::map<key, std::vector<range> > utrRangeVec;
	std::vector<std::string> transcriptVector = refGene.getKeys();
	std::map<key, std::vector<posStrandPair> > trueNegatives;
	
	//storing all UTR ranges; used to make sure a found "TN PAS" isn't potentially a PAS in another transcript
	for (auto txIter = transcriptVector.begin(); txIter != transcriptVector.end(); txIter++) {
		RefGeneProperties refGeneProp = refGene.getValueByKey(*txIter);
		utrRangeVec[refGeneProp.chr].push_back(std::make_pair(refGeneProp.txStart, refGeneProp.cdsStart));
		utrRangeVec[refGeneProp.chr].push_back(std::make_pair(refGeneProp.cdsEnd, refGeneProp.txEnd));
	}

	size_t hitCounter = 0;
	const size_t maxMatches = 30000;
	size_t inUtrCounter = 0;
	std::vector<size_t> motifFrequencies(13, 0);
	//iterating over every transcript from refGene.txt
	for (auto txIter = transcriptVector.begin(); txIter != transcriptVector.end(); txIter++) {
		RefGeneProperties refGeneProp = refGene.getValueByKey(*txIter);
		if (refGeneProp.cdsStart == refGeneProp.cdsEnd) continue;
		if (refGeneProp.chr.size() > 6 || refGeneProp.chr.empty()) continue;
		
		const std::vector<unsigned int> & txStartVector = refGeneProp.exonStarts;
		const std::vector<unsigned int> & txEndVector = refGeneProp.exonEnds;
		//iterating over all exons of a transcript
		for (size_t i = 0; i < txStartVector.size(); i++) {
			seqan::CharString temp;
			size_t idx = polar::utility::getFastaIndex(refGeneProp.chr);
			size_t start = txStartVector[i];
			size_t end = txEndVector[i];
			if (start <= refGeneProp.cdsStart) start = refGeneProp.cdsStart;
			if (end >= refGeneProp.cdsEnd) end = refGeneProp.cdsEnd;
			if (start >= end) continue;
		
			seqan::readRegion(temp, faiIndex, idx, start, end);
			if (refGeneProp.strand == "-") {
				seqan::reverseComplement(temp);
				seqan::toLower(temp);

			}
			std::string exonSequence(seqan::toCString(temp));

			auto hexamersIter = Utr3Finder::hexamers.begin();
			for (; hexamersIter != Utr3Finder::hexamers.end(); hexamersIter++) {
				auto match = exonSequence.begin();
				while ((match = std::search(match, exonSequence.end(), 
					hexamersIter->begin(), hexamersIter->end()))
					!= exonSequence.end()) {
					size_t geneticPos;
					if (refGeneProp.strand == "+") {
						geneticPos = std::distance(exonSequence.begin(), match) + start;
					} else {
						geneticPos = end - std::distance(exonSequence.begin(), match);
					}
					
					//making sure that matches aren't in the UTR of any other transcript
					bool inUtr = false;
					for (auto it = utrRangeVec[refGeneProp.chr].begin(); 
						it != utrRangeVec[refGeneProp.chr].end(); 
						it++) {
						if (geneticPos >= it->first && geneticPos <= it->second)
							inUtr = true;
							inUtrCounter++;
							break;

					}
					//searching for duplicates (overlapping transcripts/PAS positions)
					bool duplicate = false;
					for (auto it = trueNegatives[refGeneProp.chr].begin(); 
						it != trueNegatives[refGeneProp.chr].end(); 
						it++) {
						if (it->first == geneticPos) {
							duplicate = true;
							break;
						}
					}
					//check if geneticPos yields an analyzable PAS motif (correct mapping?)
					if (refGeneProp.strand == "+") {
						seqan::CharString motifAtPos;
						seqan::readRegion(motifAtPos, faiIndex, idx, geneticPos, geneticPos + 6);
						std::string stdMotifAtPos(seqan::toCString(motifAtPos));
						auto findIt = thresholdMap.find(stdMotifAtPos);
						if (findIt == thresholdMap.end()) {
							match++;
							continue;
						}
					} else {
						seqan::CharString motifAtPos;
						seqan::readRegion(motifAtPos, faiIndex, idx, geneticPos - 6, geneticPos);
						seqan::reverseComplement(motifAtPos);
						seqan::toLower(motifAtPos);
						std::string stdMotifAtPos(seqan::toCString(motifAtPos));
						auto findIt = thresholdMap.find(stdMotifAtPos);
						if (findIt == thresholdMap.end()) {
							match++;
							continue;
						}
					}

					//add to data set if found TN-PAS is neither in the UTR of a transcript 
					//nor already found in another transcript (duplicate) 
					if (! inUtr && ! duplicate) {
						trueNegatives[refGeneProp.chr].push_back(
							std::make_pair(geneticPos, refGeneProp.strand)
						);
						hitCounter++;
					}
					match++;
					if (hitCounter >= maxMatches) goto stopLoop;
				}
			}
		}
	}
	stopLoop:

	for (size_t i = 0; i < 33; i++) {
		size_t numTruePositives = 0; //found true positives/negatives
		size_t unknMotives = 0; //motives that aren't searched for or refGene.txt contains faulty mappings
		size_t totalTruePositives = 0; //total number of true positives
		size_t numTrueNegatives = 0; //true negatives
		size_t numFalsePositives = 0; //"matches" found around a true positives/negatives  (not further analyzed atm)
		size_t totalTrueNegatives = 0; //total number of true positives
		//evaluating every true positive
		std::ifstream inTP(tpFasta.string());
		std::string line;
		size_t lineCount = 0;
		while (std::getline(inTP, line)) {
			std::string strand; 
			switch (lineCount) {
			case 0:
			{	
				qi::parse(line.begin(), line.end(), 
					">" >> qi::omit[*~qi::char_('|')] >> "|" >> qi::omit[qi::uint_] >> "|" >> qi::char_, 
					strand);
				lineCount++;
			}
			case 1: 
			{
				SeqStruct ss = {
					line,
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
						break;
					} else if (resultIt->pos == 105 && resultIt->strand == "-" && strand == "-") {
						numTruePositives++;
						break;
					} else if (resultIt == u3FuzzyResVector.end() - 1) {
						//maybe write out which sequences were not correctly detected as PAS
					}
				}
				lineCount = 0;
				totalTruePositives++;
			}
			}
		}
		inTP.close()
		sensitivityVec.push_back(sensitivity);
		
		std::ifstream inTN(tnFasta.string());
		//evaluating every true negative
		while(std::getline(inTN, line)) {
			switch (lineCount) {
			case 0: 
			{
			}
			case 1: 
			{
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
				bool foundMatch = false;
				//analyzing the results
				for (auto resultIt = u3FuzzyResVector.begin(); resultIt != u3FuzzyResVector.end(); resultIt++) {
					if (resultIt->pos == 100 && resultIt->strand == "+" && strand == "+") {
						numFalsePositives++;
						foundMatch = true;
					} else if (resultIt->pos == 99 && resultIt->strand == "-" && strand == "-") {
						numFalsePositives++;
						foundMatch = true;	
					}
				}
				if (! foundMatch) {
					//std::cerr << ss.seq << std::endl;
					numTrueNegatives++;
				} else {
					//std::cerr << ss.seq << std::endl;
				}
				totalTrueNegatives++;
			}
			}
		}
		double specificity = static_cast<double>(numTrueNegatives) / totalTrueNegatives;
		specificityVec.push_back(specificity);
		Utr3FinderFuzzy tempObj(SeqStruct{
			std::string("acgt"),
			boost::none,
			boost::none,
			boost::none,
			boost::none,
			boost::none,
			boost::none,
			boost::none
			});	
		tempObj.setThresholdMap(thresholdMap);
		if (i != 0) {
			for (auto mapIter = thresholdMap.begin(); mapIter != thresholdMap.end(); mapIter++) {
				mapIter->second += 0.05;
			}
		}

	}
	std::cerr << "sensitivity,specificity" << std::endl;
	for (unsigned int i = 0; i < specificityVec.size(); i++) {
		if (i == 0) std::cerr << 0.51;
		else std::cerr << (i - 1) * 0.05;
		std::cerr << ", " 
			<< sensitivityVec[i] << "," 
			<< specificityVec[i] 
			<< std::endl;
	}
}


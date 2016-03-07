#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <chrono>
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

namespace fs = boost::filesystem;
namespace qi = boost::spirit::qi;

//writes the sequence around a putative PAS to a file (FASTA like)
bool createTNset (const fs::path & out, const seqan::FaiIndex & refGenomeIndex) {
	typedef std::string key;
	typedef std::string strand;
	typedef size_t position;
	typedef std::pair<size_t, size_t> range;
	typedef std::pair<position, strand> posStrandPair;
	
	fs::path refGeneFile = "ucsc_data/refGene.txt";
	fs::path ucscMappedTx = "ucsc_data/ucsc_txRefSeq.txt";
	if (! fs::exists(refGeneFile)) {
		std::cerr << "createTNset: " << refGeneFile.string() << " not found" << std::endl; 
		return false;
	}
	if (! fs::exists(ucscMappedTx)) {
		std::cerr << "createTNset: " << ucscMappedTx.string() << " not found" << std::endl; 
		return false;
	}

	RefGeneParser refGene(refGeneFile);	
	std::map<key, std::vector<range> > utrRangeVec;
	std::vector<std::string> transcriptVector = refGene.getKeys();
	//map that stores position and strand of a match sorted by chromosome (key)
	std::map<key, std::vector<posStrandPair> > trueNegatives;
	std::unordered_map<std::string, size_t> ucscTxRefSeq = polar::utility::getTxRefSeqAccessions(ucscMappedTx);
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

	std::unordered_map<std::string, size_t> pasCountMap = {
		{std::string("aataaa"), 0},
		{std::string("attaaa"), 0},
		{std::string("tataaa"), 0},
		{std::string("agtaaa"), 0},
		{std::string("aagaaa"), 0},
		{std::string("aatata"), 0},
		{std::string("aataca"), 0},
		{std::string("cataaa"), 0},
		{std::string("gataaa"), 0},
		{std::string("aatgaa"), 0},
		{std::string("actaaa"), 0},
		{std::string("aataga"), 0}
	};

	//storing all UTR ranges; used to make sure a found "TN PAS" isn't potentially a PAS in another transcript
	for (auto txIter = transcriptVector.begin(); txIter != transcriptVector.end(); txIter++) {
		RefGeneProperties refGeneProp = refGene.getValueByKey(*txIter);
		utrRangeVec[refGeneProp.chr].push_back(std::make_pair(refGeneProp.txStart, refGeneProp.cdsStart));
		utrRangeVec[refGeneProp.chr].push_back(std::make_pair(refGeneProp.cdsEnd, refGeneProp.txEnd));
	}

	size_t hitCounter = 0;
	const size_t maxMatches = 21600;
	const size_t maxSinglePAScount = 1800; 
	size_t inUtrCounter = 0;
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
		
			seqan::readRegion(temp, refGenomeIndex, idx, start, end);
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
					std::string stdMotifAtPos;
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
						seqan::readRegion(motifAtPos, refGenomeIndex, idx, geneticPos, geneticPos + 6);
						stdMotifAtPos = std::string(seqan::toCString(motifAtPos));
						auto findIt = thresholdMap.find(stdMotifAtPos);
						size_t numSinglePas = pasCountMap.find(stdMotifAtPos)->second;
						if (numSinglePas >= maxSinglePAScount) break;
						if (findIt == thresholdMap.end()) {
							match++;
							continue;
						}
					} else {
						seqan::CharString motifAtPos;
						seqan::readRegion(motifAtPos, refGenomeIndex, idx, geneticPos - 6, geneticPos);
						seqan::reverseComplement(motifAtPos);
						seqan::toLower(motifAtPos);
						stdMotifAtPos = std::string(seqan::toCString(motifAtPos));
						auto findIt = thresholdMap.find(stdMotifAtPos);
						size_t numSinglePas = pasCountMap.find(stdMotifAtPos)->second;
						if (numSinglePas >= maxSinglePAScount) break;
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
						pasCountMap[stdMotifAtPos]++;
					}
					match++;
					if (hitCounter >= maxMatches) goto stopLoop;
				}
			}
		}
	}
	stopLoop:
	
	std::ofstream output (out.string());
	//evaluating every true negative
	if (output.is_open()) {
		for (auto it = trueNegatives.begin(); it != trueNegatives.end(); it++) {
			//iterating over every putative PAS
			for (auto vecIt = it->second.begin(); vecIt != it->second.end(); vecIt++) {
				size_t & pos = vecIt->first;
				std::string & strand = vecIt->second;
				seqan::CharString temp;
				//mapping chromosome to id number for use in the fasta index
				unsigned idx = polar::utility::getFastaIndex(it->first);
				//copy the genomic sequence 250 bases around the genomic position of the putative PAS (500nt long)
				std::string motifSequence;
				if (strand == "+") {
					std::stringstream ss;
					seqan::readRegion(temp, refGenomeIndex, idx, pos - 250, pos + 250);
					ss << seqan::infix(temp, 250, 256);
					motifSequence = ss.str();
					seqan::erase(temp, 250,256);
				} else {
					std::stringstream ss;
					seqan::readRegion(temp, refGenomeIndex, idx, pos - 251, pos + 249);
					ss << seqan::infix(temp, 245, 251);
					motifSequence = ss.str();
					seqan::erase(temp, 245,251);
				}
				std::string seq(seqan::toCString(temp));
				unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

				std::shuffle(seq.begin(), seq.end(), std::default_random_engine(seed));
				if (strand == "+") seq.insert(250, motifSequence);
				else seq.insert(245, motifSequence);

				output << ">" << it->first << "|" << pos << "|" << strand << std::endl
					<< seq << std::endl;
			}
		}
	}
	output.close();
	return true;
}


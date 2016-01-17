#include <algorithm>
#include <climits>
#include <fstream>
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
		
	size_t numTruePositives = 0; //found true positives/negatives
	size_t unknMotives = 0; //motives that aren't searched for or refGene.txt contains faulty mappings
	size_t totalTruePositives = 0; //total number of true positives
	size_t numTrueNegatives = 0; //true negatives
	size_t numFalsePositives = 0; //"matches" found around a true positives/negatives  (not further analyzed atm)
	size_t totalTrueNegatives = 0; //total number of true positives
	std::vector<double> sensitivityVec;
	std::vector<double> specificityVec;
	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		std::cerr << "could not open index file for " << referenceGenome << std::endl;
	}
	readKnownPolyA(knownPolyA, knownPolyAvec);
	std::vector<std::pair<std::string, std::string> > utrSequenceVec;
	//creating true positive data set
	BOOST_FOREACH (KnownPolyA & knownPolyA, knownPolyAvec) {
		std::string baseSeqId;
		qi::parse(knownPolyA.id.begin(), knownPolyA.id.end(), *~qi::char_('.'), baseSeqId);
		RefGeneProperties refGeneProp = refGene.getValueByKey(baseSeqId);
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
		
		seqan::CharString utr;
		if (refGeneProp.cdsEnd < refGeneProp.exonStarts.back()) {
			if (refGeneProp.strand == "+")
				std::cerr << "Transcript: " << knownPolyA.id << "|UTR len: " << refGeneProp.txEnd - refGeneProp.cdsEnd << std::endl;
			else 
				std::cerr << "Transcript: " << knownPolyA.id << "|UTR len: " << refGeneProp.cdsStart - refGeneProp.txStart << std::endl;
			seqan::readRegion(utr, faiIndex, idx, refGeneProp.cdsEnd, refGeneProp.txEnd);
			seqan::toLower(utr);
			utrSequenceVec.push_back(std::make_pair(baseSeqId, seqan::toCString(utr)));
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
			truePositives[refGeneProp.chr].push_back(std::make_pair(geneticPos, refGeneProp.strand));
		}
	}
	//creating true negative data set
	std::map<key, std::vector<range> > utrRangeVec;
	std::vector<std::string> transcriptVector = refGene.getKeys();
	std::map<key, std::vector<posStrandPair> > trueNegatives;
	
	for (auto txIter = transcriptVector.begin(); txIter != transcriptVector.end(); txIter++) {
		RefGeneProperties refGeneProp = refGene.getValueByKey(*txIter);
		utrRangeVec[refGeneProp.chr].push_back(std::make_pair(refGeneProp.txStart, refGeneProp.cdsStart));
		utrRangeVec[refGeneProp.chr].push_back(std::make_pair(refGeneProp.cdsEnd, refGeneProp.txEnd));
	}

	size_t hitCounter = 0;
	const size_t maxMatches = 30000;
	size_t inUtrCounter = 0;
	std::vector<size_t> motifFrequencies(13, 0);
	//iterating over every available transcript from refGene.txt
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
				while ((match = std::search(match, exonSequence.end(), hexamersIter->begin(), hexamersIter->end())) 
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
	//std::cerr << "inUtrCounter: " << inUtrCounter << std::endl;
	
	std::ofstream writeTn ("trueNegatives.fa");
	if (writeTn.is_open()) {
		auto mapIter = trueNegatives.begin();
		for (; mapIter != trueNegatives.end(); mapIter++) {
			for (auto vecIter = mapIter->second.begin(); vecIter != mapIter->second.end(); vecIter++) {
				seqan::CharString temp;
				size_t idx = polar::utility::getFastaIndex(mapIter->first);
				seqan::readRegion(temp, faiIndex, idx, vecIter->first - 100, vecIter->first + 100);
				writeTn << ">" << mapIter->first << "|" << vecIter->first << "|" << vecIter->second << std::endl;
				writeTn << temp << std::endl; 
			}
		}
	}
	writeTn.close();

	std::ofstream writeTp ("truePositivesUTR.fa");
	if (writeTp.is_open()) {
		auto vecIter = utrSequenceVec.begin();
		for (; vecIter != utrSequenceVec.end(); vecIter++) {
			writeTp << ">UTR of: " << vecIter->first << std::endl;
			writeTp << vecIter->second << std::endl;
		}
	}
	writeTp.close();
	utrSequenceVec.clear();

	for (size_t i = 0; i < 32; i++) {
		numTruePositives = 0; //found true positives/negatives
		unknMotives = 0; //motives that aren't searched for or refGene.txt contains faulty mappings
		totalTruePositives = 0; //total number of true positives
		numTrueNegatives = 0; //true negatives
		numFalsePositives = 0; //"matches" found around true positives/negatives  (not further analyzed atm)
		totalTrueNegatives = 0; //total number of true positives
		
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
					} else if (resultIt->pos == 105 && resultIt->strand == "-" && strand == "-") {
						numTruePositives++;
					}
				}
				std::string motif = ss.seq.substr(100, 6);
				if (strand == "-") {
					std::string temp;
					std::transform(motif.rbegin(), motif.rend(), 
						std::back_inserter(temp), 
						polar::utility::complement);
					motif = temp;
				}
				
				//std::cerr << resultIt->truthValue << "  " << motif << std::endl;
				if (Utr3FinderFuzzy::dseLocMap.find(motif) == Utr3FinderFuzzy::dseLocMap.end()) {
					unknMotives++;			
				}
				totalTruePositives++;
			}
		}
		//std::cerr << "Threshold: " << thresholdMap.find("aataaa")->second << std::endl;
		//std::cerr << "-------Sensitivity test-------" << std::endl;
		//std::cerr << "correct predictions (found true positives): " << numTruePositives << std::endl;
		//std::cerr << "incorrect predictions (false negatives): " << totalTruePositives - unknMotives - numTruePositives << std::endl;
		//std::cerr << "total sequences analyzed (total true positives): " << totalTruePositives << std::endl; 
		//std::cerr << "unknown motives: " << unknMotives << std::endl;
		double sensitivity = static_cast<double>(numTruePositives) / (totalTruePositives - unknMotives);
		sensitivityVec.push_back(sensitivity);

		//evaluating every true negative
		for (auto it = trueNegatives.begin(); it != trueNegatives.end(); it++) {
			//iterating over every putative PAS
			for (auto vecIt = it->second.begin(); vecIt != it->second.end(); vecIt++) {
				size_t & pos = vecIt->first;
				std::string & strand = vecIt->second;
				seqan::CharString temp;
				//mapping chromosome to id number for use in the fasta index
				unsigned idx = polar::utility::getFastaIndex(it->first);
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
					//std::cerr << ">" << pos << "|" << strand << "|" << it->first << std::endl;
					//std::cerr << ss.seq << std::endl;
				}
				totalTrueNegatives++;
			}
		}
		//for (auto it = motifFrequencies.begin(); it != motifFrequencies.end(); it++) {
		//	std::cerr << *it << " ";
		//}
//		std::cerr << "-------Specifity test-------" << std::endl;
//		std::cerr << "correct predictions (found true negatives): " << numTrueNegatives << std::endl;
//		std::cerr << "incorrect predictions (false positives): " << numFalsePositives << std::endl;
//		std::cerr << "total sequences analyzed (total true positives): " << totalTrueNegatives << std::endl; 
		//std::cerr << "incorrect predictions: " << numFalsePositives << std::endl;
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
		for (auto mapIter = thresholdMap.begin(); mapIter != thresholdMap.end(); mapIter++) {
			mapIter->second += 0.05;
		}


	}
	std::cerr << "totalTruePositives: " << totalTruePositives << std::endl;
	std::cerr << "sensitivity,specificity" << std::endl;
	for (unsigned int i = 0; i < specificityVec.size(); i++) {
		std::cerr << sensitivityVec[i] << "," << specificityVec[i] << std::endl;
	}
}


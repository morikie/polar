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
#include "readKnownPolyA.hpp"
#include "createTPset.hpp"

namespace fs = boost::filesystem;
namespace qi = boost::spirit::qi;


bool createTPset (const fs::path & out, const seqan::FaiIndex & refGenomeIndex) {
	typedef std::string key;
	typedef std::string strand;
	typedef size_t position;
	typedef std::pair<position, strand> posStrandPair;
	//map that stores position and strand of a match sorted by chromosome (key)
	std::map<key, std::vector<posStrandPair> > truePositives;
	fs::path refGeneFile = "ucsc_data/refGene.txt";
	fs::path knownPolyA = "../perf_testing/knownPolyAtranscript.txt";
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path ucscMappedTx = "ucsc_data/ucsc_txRefSeq.txt";
	std::unordered_map<std::string, size_t> ucscTxRefSeq = polar::utility::getTxRefSeqAccessions(ucscMappedTx);
	std::vector<KnownPolyA> knownPolyAvec;
	RefGeneParser refGene(refGeneFile);	
	
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
		
		//calculate size of the transcript by using exon start/end values
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
				seqan::readRegion(motifAtPos, refGenomeIndex, idx, geneticPos, geneticPos + 6);
				std::string stdMotifAtPos(seqan::toCString(motifAtPos));
				auto findIt = thresholdMap.find(stdMotifAtPos);
				if (findIt == thresholdMap.end()) {
					continue;
				}
			} else {
				seqan::CharString motifAtPos;
				seqan::readRegion(motifAtPos, refGenomeIndex, idx, geneticPos - 6, geneticPos);
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
	std::ofstream output(out.string());
	if (output.is_open()) {
		for (auto it = truePositives.begin(); it != truePositives.end(); it++) {
			for (auto vecIt = it->second.begin(); vecIt != it->second.end(); vecIt++) {
				size_t & pos = vecIt->first;
				std::string & strand = vecIt->second;
				seqan::CharString temp;
				size_t idx = polar::utility::getFastaIndex(it->first);
				//copy the genomic sequence 100 bases around the genomic position of the PAS (200nt long)
				if (strand == "+") {
					seqan::readRegion(temp, refGenomeIndex, idx, pos - 250, pos + 250);
				} else {
					seqan::readRegion(temp, refGenomeIndex, idx, pos - 251, pos + 249);
				}
				std::string seq(seqan::toCString(temp));
				
				output << ">" << it->first << "|" << pos << "|" << strand << std::endl;
				output << seq << std::endl;
				
			}
		}
	}
	output.close();
	return true;
}

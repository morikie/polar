#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/spirit/include/qi.hpp>
#include <seqan/seq_io.h>
#include "../src/hgvsParser.hpp"
#include "../src/polarUtility.hpp"
#include "../src/utr3Finder.hpp"
#include "../src/utr3FinderFuzzy.hpp"
#include "../src/utr3FinderNaive.hpp"
#include "createTNset.hpp"
#include "createTPset.hpp"
#include "perf_fuzzy.hpp"

namespace fs = boost::filesystem;
namespace qi = boost::spirit::qi;

//structure needed for parsing the data sets
struct PasPosition {
	std::string chr;
	size_t pos;
	std::string strand;

	void reset() {
		this->chr = std::string();
		this->pos = 0;
		this->strand = std::string();
	}
};

//used by boost::spirit to parse into a PasPosition object
BOOST_FUSION_ADAPT_STRUCT (
	PasPosition,
	(std::string, chr)
	(size_t, pos)
	(std::string, strand)
)

int main (int argc, char * argv[]) {
	//map that stores position and strand of a match sorted by chromosome (key)
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	fs::path ucscMappedTx = "ucsc_data/ucsc_txRefSeq.txt";
	fs::path tpFasta = "../perf_testing/tpSet.fa";
	fs::path tnFasta = "../perf_testing/tnSet.fa";
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
	std::vector<double> sensitivityVec;
	std::vector<double> specificityVec;
	size_t totalTrueNegatives = 0; //total number of true negatives
	size_t totalTruePositives = 0; //total number of true positives
	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		std::cerr << "could not open index file for " << referenceGenome << std::endl;
	}
	if (! fs::exists(tpFasta)) {
		createTPset(tpFasta, faiIndex);
	} 
	if (! fs::exists(tnFasta)) {
		createTNset(tnFasta, faiIndex);
	}
	
	//data structures used to store the TP/TN data sets; stores chromosome, position, strand (PasPosition) and 
	//a 2x 250nt long sequence around the position (500 in total)
	std::vector<std::pair<PasPosition, SeqStruct> > pasPosAndSeqTP;
	std::vector<std::pair<PasPosition, SeqStruct> > pasPosAndSeqTN;	
	PasPosition pasPos;
	std::ifstream inTP(tpFasta.string());
	std::string line;
	size_t lineCount = 0;
	while (std::getline(inTP, line)) {
		switch (lineCount) {
		case 0:
		{	
			qi::parse(line.begin(), line.end(), 
				">" >> +~qi::char_('|') >> "|" >> qi::uint_ >> "|" >> +qi::char_, 
				pasPos);
			lineCount++;
			break;
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
			//storing position, chromosome and strand of the TP and the sequence around the TP
			pasPosAndSeqTP.push_back(std::make_pair(pasPos, ss));

			pasPos.reset();
			lineCount = 0;
			totalTruePositives++;
			break;
		}
		}
	}
	inTP.close();
	
	//reading in the TN data set
	std::ifstream inTN(tnFasta.string());
	lineCount = 0;
	line.clear();
	pasPos.reset();
	while(std::getline(inTN, line)) {
		switch (lineCount) {
		case 0: 
		{
			qi::parse(line.begin(), line.end(), 
				">" >> +~qi::char_('|') >> "|" >> qi::uint_ >> "|" >> +qi::char_, 
				pasPos);
			lineCount++;
			break;
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
			//storing position, chromosome and strand of the TN and the sequence around the TN
			pasPosAndSeqTN.push_back(std::make_pair(pasPos, ss));
			
			totalTrueNegatives++;
			pasPos.reset();
			lineCount = 0;
			break;
		}
		}
	}
	inTN.close();

	//running the prediction 33 times with different thresholds (0 to 1.6 in 0.05 steps and one more for the default thresholds)
	for (size_t i = 0; i < 33; i++) {
		size_t numTruePositives = 0; //found true positives/negatives
		size_t numTrueNegatives = 0; //true negatives
		size_t numFalsePositives = 0; //"matches" found around a true positive/negative
		double sensitivity = 0.0;
		double specificity = 0.0;
		
		//evaluating every true positive
		for (auto vecIt = pasPosAndSeqTP.begin(); vecIt != pasPosAndSeqTP.end(); vecIt++) { 
			Utr3FinderFuzzy u3Fuzzy(vecIt->second);		
			std::vector<Utr3Finder::Utr3FinderResult> u3FuzzyResVector = u3Fuzzy.getPolyaMotifPos();
			std::string & strand = vecIt->first.strand;

			//analyzing the results
			for (auto resultIt = u3FuzzyResVector.begin(); resultIt != u3FuzzyResVector.end(); resultIt++) {
				//std::cerr << resultIt->pos << "(" << resultIt->strand << "|" << strand << "), ";
				if (resultIt->pos == 250 && resultIt->strand == strand) {
					numTruePositives++;
					break;
				} else {
					//maybe write out which sequences were not correctly detected as PAS
				}
			}
			//std::cerr << std::endl;
		}
		sensitivity = static_cast<double>(numTruePositives) / totalTruePositives;
		sensitivityVec.push_back(sensitivity);
		
		//evaluating every true negative
		for (auto vecIt = pasPosAndSeqTN.begin(); vecIt != pasPosAndSeqTN.end(); vecIt++) {
			Utr3FinderFuzzy u3Fuzzy(vecIt->second);
			std::vector<Utr3Finder::Utr3FinderResult> u3FuzzyResVector = u3Fuzzy.getPolyaMotifPos();
			std::string & strand = vecIt->first.strand;
			bool foundMatch = false;
			
			//analyzing the results
			for (auto resultIt = u3FuzzyResVector.begin(); resultIt != u3FuzzyResVector.end(); resultIt++) {
				//std::cerr << resultIt->pos << "(" << resultIt->strand << "|" << strand << "), ";
				if (resultIt->pos == 250 && resultIt->strand == strand) {
					numFalsePositives++;
					foundMatch = true;
				}
			}
			//std::cerr << std::endl;
			if (! foundMatch) {
				//std::cerr << ss.seq << std::endl;
				numTrueNegatives++;
			} else {
				//std::cerr << ss.seq << std::endl;
			}
		}

		
		specificity = static_cast<double>(numTrueNegatives) / totalTrueNegatives;
		specificityVec.push_back(specificity);

		//temporary object to alter the static threshold map
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
	
	std::cerr << "threshold,sensitivity,specificity" << std::endl;
	for (unsigned int i = 0; i < specificityVec.size(); i++) {
		if (i == 0) std::cerr << 0.51;
		else std::cerr << (i - 1) * 0.05;
		std::cerr << "," 
			<< sensitivityVec[i] << "," 
			<< specificityVec[i] 
			<< std::endl;
	}
}


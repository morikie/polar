#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <map>
#include <omp.h>
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
	typedef std::unordered_map<std::string, size_t> map;
	//map that stores position and strand of a match sorted by chromosome (key)
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	fs::path ucscMappedTx = "ucsc_data/ucsc_txRefSeq.txt";
	fs::path canonicalTpFasta = "../scripts/polyadq_test/canonicalTpSet.fa";
	fs::path canonicalTnFasta = "../scripts/polyadq_test/canonicalTnSet.fa";
	fs::path positiveFasta = "../perf_testing/positiveSet.fa";
	fs::path negativeFasta = "../perf_testing/negativeSet.fa";
	fs::path tpFasta = "TPdataSet.fa";
	fs::path fnFasta = "FNdataSet.fa";
	fs::path tnFasta = "TNdataSet.fa";
	fs::path fpFasta = "FPdataSet.fa";
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
	std::unordered_map<std::string, size_t> pasDistributionTN = {
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
	std::unordered_map<std::string, size_t> pasDistributionTP = {
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
	std::unordered_map<std::string, std::ofstream*> pasFileMap = {
		{std::string("aataaa"), new std::ofstream("aataaa.fa")},
		{std::string("attaaa"), new std::ofstream("attaaa.fa")},
		{std::string("tataaa"), new std::ofstream("tataaa.fa")},
		{std::string("agtaaa"), new std::ofstream("agtaaa.fa")},
		{std::string("aagaaa"), new std::ofstream("aagaaa.fa")},
		{std::string("aatata"), new std::ofstream("aatata.fa")},
		{std::string("aataca"), new std::ofstream("aataca.fa")},
		{std::string("cataaa"), new std::ofstream("cataaa.fa")},
		{std::string("gataaa"), new std::ofstream("gataaa.fa")},
		{std::string("aatgaa"), new std::ofstream("aatgaa.fa")},
		{std::string("actaaa"), new std::ofstream("actaaa.fa")},
		{std::string("aataga"), new std::ofstream("aataga.fa")}
	};
	std::vector<double> sensitivityVec;
	std::vector<double> specificityVec;
	size_t totalTrueNegatives = 0; //total number of true negatives
	size_t totalTruePositives = 0; //total number of true positives
	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		std::cerr << "could not open index file for " << referenceGenome << std::endl;
	}
	if (! fs::exists(positiveFasta)) {
		createTPset(positiveFasta, faiIndex);
	} 
	if (! fs::exists(negativeFasta)) {
		createTNset(negativeFasta, faiIndex);
	}
	
	//data structures used to store the TP/TN data sets; stores chromosome, position, strand (PasPosition) and 
	//a 2x 250nt long sequence around the position (500 in total)
	std::vector<std::pair<PasPosition, SeqStruct> > pasPosAndSeqTP;
	std::vector<std::pair<PasPosition, SeqStruct> > pasPosAndSeqTN;	
	PasPosition pasPos;
	std::ifstream inTP(positiveFasta.string());
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
			
			std::string pas;
			if (pasPos.strand == "-") {
				std::transform(line.rbegin() + 249, line.rbegin() + 255, 
					std::back_inserter(pas), 
					polar::utility::complement);
			} else {
				pas = std::string(line.begin() + 250, line.begin() + 256);
			}
			if (pasDistributionTP.find(pas) != pasDistributionTP.end()) {
				pasDistributionTP[pas]++;
			} else {
				std::cerr << pas << " was not found in TP set!" << std::endl;		
			}
			*(pasFileMap.find(pas)->second) << ">" <<  pasPos.chr << "|" 
				<< pasPos.pos << "|" 
				<< pasPos.strand << std::endl
				<< line << std::endl;
			pasPos.reset();
			lineCount = 0;
			totalTruePositives++;
			break;
		}
		}
	}
	inTP.close();
	
	//reading in the TN data set
	std::ifstream inTN(negativeFasta.string());
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
			
			std::string pas;
			if (pasPos.strand == "-") {
				std::transform(line.rbegin() + 249, line.rbegin() + 255, 
					std::back_inserter(pas), 
					polar::utility::complement);
			} else {
				pas = std::string(line.begin() + 250, line.begin() + 256);
			}
			
			if (pasDistributionTN.find(pas) != pasDistributionTN.end()) {
				pasDistributionTN[pas]++;
			} else {
				std::cerr << pas << " was not found in TN set!" << std::endl;		
			}

			totalTrueNegatives++;
			pasPos.reset();
			lineCount = 0;
			break;
		}
		}
	}
	std::cerr << "TP pas distribution" << std::endl;
	BOOST_FOREACH(map::value_type & pair, pasDistributionTP) {
		std::cerr << pair.first << ": " << pair.second << std::endl;
	}

	std::cerr << "TN pas distribution" << std::endl;
	BOOST_FOREACH(map::value_type & pair, pasDistributionTN) {
		std::cerr << pair.first << ": " << pair.second << std::endl;
	}
	inTN.close();
	bool writeDataSets = false;
	std::ofstream outTP(tpFasta.string(), std::ofstream::out);
	std::ofstream outFN(fnFasta.string(), std::ofstream::out);
	std::ofstream outTN(tnFasta.string(), std::ofstream::out);
	std::ofstream outFP(fpFasta.string(), std::ofstream::out);
	//running the prediction 33 times with different thresholds (0 to 1.6 in 0.05 steps and one more for the default thresholds)
	double step = 0.05;
	for (size_t i = 0; i < 45; i++) {
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
			bool foundMatch = false;

			//analyzing the results
			for (auto resultIt = u3FuzzyResVector.begin(); resultIt != u3FuzzyResVector.end(); resultIt++) {
				//std::cerr << resultIt->pos << "(" << resultIt->strand << "|" << strand << "), ";
				if (resultIt->pos == 250 && resultIt->strand == strand) {
					//TruePositives
					if (writeDataSets) {
						outTP << ">" <<  vecIt->first.chr << "|" 
							<< vecIt->first.pos << "|" 
							<< strand << std::endl
							<< u3Fuzzy.getSequence() << std::endl;
					}
					foundMatch = true;
					numTruePositives++;
					break;
				}
			}
			//std::cerr << std::endl;
			if (! foundMatch) {
				//FalseNegatives
				if (writeDataSets) {
					outFN << ">" <<  vecIt->first.chr << "|" 
						<< vecIt->first.pos << "|" 
						<< strand << std::endl
						<< u3Fuzzy.getSequence() << std::endl;
				}
			}
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
					//FalsePositives
					if (writeDataSets) {
						outFP << ">" <<  vecIt->first.chr << "|" 
							<< vecIt->first.pos << "|" 
							<< strand << std::endl
							<< u3Fuzzy.getSequence() << std::endl;
					}
					numFalsePositives++;
					foundMatch = true;
					break;
				}
			}
			//std::cerr << std::endl;
			if (! foundMatch) {
				//TrueNegatives
				if (writeDataSets) {
					outTN << ">" <<  vecIt->first.chr << "|" 
						<< vecIt->first.pos << "|" 
						<< strand << std::endl
						<< u3Fuzzy.getSequence() << std::endl;
				}
				numTrueNegatives++;
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
				mapIter->second += step;
			}
		}
	}
	outTN.close();
	outTP.close();
	outFN.close();
	outFP.close();
	std::cerr << "threshold,sensitivity,specificity" << std::endl;
	for (unsigned int i = 0; i < specificityVec.size(); i++) {
		if (i == 0) std::cerr << "default";
		else std::cerr << (i - 1) * step;
		std::cerr << "," 
			<< sensitivityVec[i] << "," 
			<< specificityVec[i] 
			<< std::endl;
	}
}


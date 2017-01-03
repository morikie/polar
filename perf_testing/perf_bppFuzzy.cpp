#include <iostream>
#include <string>
#include <boost/filesystem/path.hpp>
#include <boost/foreach.hpp>
#include <boost/spirit/include/qi.hpp>
#include "../src/bppPredictFuzzy.hpp"
#include "../src/polarUtility.hpp"
#include "../src/refGeneParser.hpp"
#include "../src/utr3FinderFuzzy.hpp"
#include "perf_bppFuzzy.hpp"


namespace fs = boost::filesystem;
namespace qi = boost::spirit::qi;

//structure needed for parsing the data sets
struct PasPosition {
	std::string seqId;
	std::string chr;
	size_t pos;
	std::string strand;

	void reset() {
		this->chr = std::string();
		this->pos = 0;
		this->strand = std::string();
		this->seqId = std::string();
	}
};

//used by boost::spirit to parse into a PasPosition object
BOOST_FUSION_ADAPT_STRUCT (
	PasPosition,
	(std::string, seqId)
	(std::string, chr)
	(size_t, pos)
	(std::string, strand)
)


/**
 * Parsing FASTA-like formats
 */
void getPositiveDataSet(std::vector<PasPosition> & returnVec, const fs::path & positiveFasta) {
	std::ifstream in(positiveFasta.string());
	std::string line;
	size_t lineCount = 0;
	PasPosition pasPos;

	while (std::getline(in, line)) {
		switch (lineCount) {
		case 0:
		{	
			qi::parse(line.begin(), line.end(), 
				">" >> +~qi::char_('.') >> qi::omit[*~qi::char_("|")] >> "|" 
				>> +~qi::char_('|') >> "|" 
				>> qi::uint_ >> "|" 
				>> +qi::char_, 
				pasPos);
			returnVec.push_back(pasPos);
			pasPos.reset();
			lineCount++;
			break;
		}
		case 1: 
		{
			lineCount = 0;
			break;
		}
		}
	}
	in.close();
}


/**
 * Parse all transcripts and fill up a map containing the UTR sequence of the transcript.
 */
void getUtrByTranscript(std::unordered_map<std::string, std::string> & returnMap, fs::path & path) {
	std::ifstream in(path.string());
	std::string line;
	size_t lineCount = 0;
	std::string tempString;

	while (std::getline(in, line)) {
		switch (lineCount) {
		case 0:
		{	
			qi::parse(line.begin(), line.end(), 
				">" >> +~qi::char_('.') >> "." >> qi::omit[*qi::char_], 
				tempString);
			++lineCount;
			break;
		}
		case 1: 
		{
			if (returnMap.find(tempString) == returnMap.end()) {
				//std::cerr << __LINE__ << std::endl;
				returnMap.insert(std::pair<std::string, std::string>(tempString, line));

			}
			tempString.clear();
			lineCount = 0;
			break;
		}
		}
	}
	in.close();
}


bool testCase(RefGeneParser & parser) {
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		std::cerr << "could not open index file for " << referenceGenome << std::endl;
	}
	std::string testId = "NM_004081";
	std::string testUtr = "taaattccgttgttactcaagatgactgcttcaagggtaaaagagtgcatcgctttagaagaagtttggcagtatttaaatctgttggatcctctcagctatctagtttcatgggaagttgctggttttgaatattaagctaaaagttttccactattacagaaattctgaattttggtaaatcacactgaaactttctgtataacttgtattattagactctctagttttatcttaacactgaaactgttcttcattagatgtttatttagaacctggttctgtgtttaatatatagtttaaagtaacaaataatcgagactgaaagaatgttaagatttatctgcaaggatttttaaaaaattgaaacttgcattttaagtgtttaaaagcaaatactgactttcaaaaaagtttttaaaacctgatttgaaagctaacaattttgatagtctgaacacaagcatttcacttctccaagaagtacctgtgaacagtacaatatttcagtattgagctttgcatttatgatttatctagaaatttacctcaaaagcagaatttttaaaactgcatttttaatcagtggaactcaatgtatagttagctttattgaagtcttatccaaacccagtaaaacagattctaagcaaacagtccaatcagtgagtcataatgtttattcaaagtattttatcttttatctagaatccacatatgtatgtccaatttgattgggatagtagttaggataactaaaattctgggcctaattttttaaagaatccaagacaaactaaactttactgggtatataaccttctcaatgagttaccattcttttttataaaaaaaattgttccttgaaatgctaaacttaatggctgtatgtgaaatttgcaaaatactggtattaaagaacgctgcagcttttttatgtcactcaaaggttaatcggagtatctgaaaggaattgtttttataaaaacattgaagtattagttacttgctataaatagatttttatttttgttttttagcctgttatatttccttctgtaaaataaaatatgtccagaagaggcatgttgtttctagattaggtagtgtcctcattttatattgtgaccacacagctagagcaccagagcccttttgctatactcacagtcttgttttcccagcctcttttactagtctttcaggaggtttgctcttagaactggtgatgtaaagaatggaagtagctgtatgagcagttcaaaggccaagccgtggaatggtagcaatgggatataatacctttctaagggaaacatttgtatcagtatcatttgatctgccatggacatgtgtttaaagtggctttctggcccttctttcaatggcttcttccctaaaacgtggagactctaagttaatgtcgttactatgggccatattactaatgcccactggggtctatgatttctcaaaattttcattcggaatccgaaggatacagtctttaaactttagaattcccaagaaggctttattacacctcagaaattgaaagcaccatgactttgtccattaaaaaattatccatagtttttttagtgcttttaacattccgacatacatcattctgtgattaaatctccagatttctgtaaatgatacctacattctaaagagttaattctaattattccgatatgaccttaaggaaaagtaaaggaataaatttttgtctttgttgaagtatttaatagagtaaggtaaagaagatattaagtccctttcaaaatggaaaattaattctaaactgagaaaaatgttcctactacctattgctgatactgtctttgcataaatgaataaaaataaactttttttcttcaaatgtg";
	std::string testOffset = "tttttggctttccgatgtaataatgtaaaatggtggggagttgcgtgggaactgtgtaacaaggtttaaattcgtataacaagctttagattcttaaaatgcagaagtataaagttcagtatactaatctgtctgagttagcccataaaagcaaatgtaggtacaaagataagtttaagaggtgcatcaacagcagtgcag";
	auto txPro = parser.getValueByKey(testId);
	std::string extractedUtr = polar::utility::getUtrSequence(txPro, faiIndex);
	
	size_t testLengthFromGenome = polar::utility::getUtrLength(txPro);
	//std::cerr << "utrLenFromGenome: " << testLengthFromGenome << " | utrLen NM_004081: " << testUtr.size() << std::endl;
	if (extractedUtr != testUtr) {
		std::cerr << "extractedUtr:" << std::endl << extractedUtr << std::endl << testId << std::endl << testUtr << std::endl;
		return false;
	}
	size_t zero = 25345239;
	size_t two = 25345237;
	size_t shouldBeZero = polar::utility::mapGenomePosToTxPos(txPro, zero);
	size_t shouldBeTwo = polar::utility::mapGenomePosToTxPos(txPro, two);
	if (shouldBeZero != 0 || shouldBeTwo != 2) {
		std::cerr << "shouldBeZero: " << shouldBeZero << std::endl;
		std::cerr << "shouldBeTwo: " << shouldBeTwo << std::endl;
		//return false;
	}
	BppPredictFuzzy tempObj = BppPredictFuzzy(testId);
	std::string tempObjUtr = tempObj.getUtrSeq();
	std::string tempObjOffsetSeq = tempObj.getOffsetSeq();
	if (testUtr != tempObjUtr) return false;
	if (tempObjOffsetSeq != testOffset) {
		std::cerr << testOffset << std::endl 
			<< tempObjOffsetSeq << std::endl;
		return false;
	}
	std::string testId2 = "NM_018836";
	std::string testUtr2 = "ctggccgaagtcttttttacctcctgggggcagggcagacgccgtgtgtctgtttcacggattccgttggtgaacctgtaaaaacaaaacaaacaaaacaaaacaaaaaagacaaaacctaaaactgagctatctaagggggagggtccccgcacctaccacttctgtttgccggtgggaaactcacagagcaggacgctctaggccaaatctatttttgtaaaaatgctcatgcctatgggtgactgccttctcccagagttttctttggagaacagaaagaagaaaggaaagaaaggaaccagaggcagagagacgaggatacccagcgaaagggacgggaggaagcatccgaaacctaggattcgtcctacgattctgaacctgtgccaataataccattatgtgccatgtactgacccgaaaggctcggccgcagagccggggcccagcgaatcacgcagagaaatcttacagaaaacaggggtgggaatctcttccgatagagtcgctatttctggttaatatacatatataaatatataaatacaaacacacacacacactttttttgtactgtagcaatttttgaagatcttaaatgttcctttttaaaaaaaagaattgtgttataggttacaaaatctgatttatttaacatgcttagtatgagcagaataaaccagtgttttctactttggcaactcacgtcacacacatattacacacatgtgcgcattacacacacacaatacacatacatgcatatagacgcatctattggaaatgcagttccacaggtgagcatgttctttctggtgacctggtattccatcaccattcaccccaggggacagcctcgaccgagacaaggaggcccttaaatgacagcctgcatttgctagacggttggtgagtggcatcaaatgtgtgacttactatcttgggccagaactaagaatgccaaggttttatatatgtgtgtgtatatatatatatatatatatatatatatatatatatatgtttgtgtgtgtatatatatatatatatatatatgtttgtgtgtgtatatatatgtttgtgtatatatatacacatatgcatacatatgatttttttttttcatttaagtgttggaagatgctacctaacagccacgttcacatttacgtagctggttgcttacaaacgggcctgagcccctggttgggtgggtggtggattcttggacgtgtgtgtcatacaagcatagactggattaaagaagttttccagttccaaaaattaaaggaatatatcctta";
	std::string testOffset2 = "tgatgtgtgtgtgtaatatcagggcagaacttagacatacgtgaagggccccggttggtttgaaaacgaaaaatagtcattctgtgtgcaaaccacaaggctgccccagtcaggcagcgccctgacctggcctgtgctgcattgccttcccttgcgcaggtgggcaggtgtggcccgctttttctagggcccaagggtgac";
	txPro = parser.getValueByKey(testId2);
	extractedUtr = polar::utility::getUtrSequence(txPro, faiIndex);
	
	zero = 4715104;
	two = 4715106;
	shouldBeZero = polar::utility::mapGenomePosToTxPos(txPro, zero);
	shouldBeTwo = polar::utility::mapGenomePosToTxPos(txPro, two);

	if (shouldBeZero != 0 || shouldBeTwo != 2) {
		std::cerr << "shouldBeZero: " << shouldBeZero << std::endl;
		std::cerr << "shouldBeTwo: " << shouldBeTwo << std::endl;
		return false;
	}
	BppPredictFuzzy tempObj2 = BppPredictFuzzy(testId2);
	std::string tempObj2Utr = tempObj2.getUtrSeq();
	std::string tempObj2OffsetSeq = tempObj2.getOffsetSeq();

	if (testUtr2 != tempObj2Utr) {
		std::cerr << testUtr2 << std::endl 
			<< tempObj2Utr << std::endl;
		return false;
	}
	if (tempObj2OffsetSeq != testOffset2) {
		std::cerr << testOffset2 << std::endl 
			<< tempObj2OffsetSeq << std::endl;
		return false;
	}
	
	testLengthFromGenome = polar::utility::getUtrLength(txPro);
	//std::cerr << "utrLenFromGenome: " << testLengthFromGenome << " | utrLen " << testId2 << ": " << testUtr2.size() << std::endl;
	if (extractedUtr != testUtr2) {
		std::cerr << "extractedUtr:" << std::endl << extractedUtr << std::endl << testId2 << std::endl << testUtr2 << std::endl;
		return false;
	}
	return true;

}


int main (int argc, char * argv[]) {
	//TODO: scripts/positiveSet.fa und scripts/negativeSet.fa einlesen, dann bearbeiten mit BppPredictFuzzy und Sens/Spec berechnen
	//am besten lange Sequenzen rausfiltern
	fs::path bppOut = "../perf_testing/bppOutput.txt";
	fs::path positiveDataset = "../perf_testing/positiveSet.fa";
	fs::path utrSequences = "../scripts/rna_UTRs_copy.fa";
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	fs::path refGenPath = "ucsc_data/refGene.txt";
	std::vector<PasPosition> pasVector;
	std::unordered_map<std::string, std::string> utrSeq;	
	getPositiveDataSet(pasVector, positiveDataset);
	getUtrByTranscript(utrSeq, utrSequences);
	size_t utrSeqLength = utrSeq.size();
	double tvThreshold = 1.0;
	typedef size_t truePos;
	std::vector<std::pair<BppPredictFuzzy, truePos> > resVector;
	RefGeneParser refGen = RefGeneParser(refGenPath);		

	if (! testCase(refGen)) return 1;		
	
	
	std::unordered_map<Utr3FinderFuzzy::motifSequence, double> thresholdMap = {
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
	Utr3FinderFuzzy::setThresholdMap(thresholdMap);
	size_t ignoredPas = 0;	
	size_t loopCounter = 0;
	size_t evaluatedPas = 0;
	std::set<std::string> uniques;
	std::ofstream bppOutStream(bppOut.string());
	for (auto & pas : pasVector) {
		auto txProp = refGen.getValueByKey(pas.seqId);
		size_t utrLength = polar::utility::getUtrLength(txProp);
		size_t txPos = polar::utility::mapGenomePosToTxPos(txProp, pas.pos);
		size_t txLen = polar::utility::getTxLength(txProp);
		if (txPos == UINT_MAX) {
			std::cerr << "warning: couldn't map positions for " << pas.seqId << std::endl;
		}
		if (static_cast<int>(txPos) - (static_cast<int>(txLen) - static_cast<int>(utrLength)) < 0 ) {
			std::cerr << "uint overflow: " << pas.seqId << " | txPos: " << txPos
				<< " | txLen: " << txLen 
				<< " | utrLength: " << utrLength << std::endl;
		}	
		if (argc > 1 &&std::string(argv[1]) == "verbose") {
			std::cerr << pas.seqId << " | " << "utrLength: " << utrLength  << " | ";
		}
		
		if (utrLength < 5000 && txPos != UINT_MAX) {
			BppPredictFuzzy tempObj = BppPredictFuzzy(pas.seqId);
			resVector.push_back(std::make_pair(tempObj, txPos - (txLen - utrLength)));
			evaluatedPas++;
			if (uniques.insert(tempObj.getTxId()).second) {
				std::vector<double> bppVector = tempObj.getMaxBppVector();
				bppOutStream << ">" << tempObj.getTxId() << "\n";
				for(double & d : bppVector) {
					bppOutStream << d << " ";
				}
				bppOutStream << "\n";
			}
		} else {
			++ignoredPas;
		}
		if (argc > 1 && std::string(argv[1]) == "verbose") {
			double process = static_cast<double>(loopCounter) / pasVector.size() * 100;
			std::cout << "Process: " <<  std::setprecision(2) << std::setw(6) << std::fixed << process << "%" << std::endl;
				//<< "\r" << std::flush;
		}
		++loopCounter;
	}
	
	std::cout << std::endl << "ignoredPas: " << ignoredPas << std::endl;
	std::cout << "evaluatedPas: " << evaluatedPas << std::endl;
	
	size_t total = resVector.size();
	double threshold = 1.0;
	double stepSize = 0.05;
	//sensitivity/specificity evaluation
	for (;threshold <= 5.5; threshold += stepSize) {
		size_t notFound = 0;
		std::ofstream sensitivityOut("sensitivityOut.csv");
		for (auto & resObject : resVector) {
			auto posObjVector = resObject.first.getResults(threshold);
			bool found = false;
			for (auto & posObj : posObjVector) {
				if (posObj.pos == resObject.second) {
					sensitivityOut << resObject.first.getTxId() << "|" 
						<< posObj.pos << ", "
						<< posObj.truthValue << std::endl;
					found = true;
					break;
				}
			}
			if (! found) {
				//std::cerr << resObject.first.getTxId() <<  
				//	" true position: " << resObject.second << std::endl;
				notFound++;
			}
		}	
		sensitivityOut.close();
		double sensitivity = (static_cast<double>(total) - notFound) / total;
		std::cerr << "Sensitivity at threshold: " << threshold << " | " << std::fixed << std::setprecision(2) << sensitivity << std::endl;
		//std::cerr << "total: " << total << std::endl;
		//std::cerr << "not found: " << notFound << std::endl;
	}
	return EXIT_SUCCESS;
}


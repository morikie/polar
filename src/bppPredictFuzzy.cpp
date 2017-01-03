#include <string>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/spirit/include/qi.hpp>
#include <seqan/seq_io.h>
#include "polarUtility.hpp"
#include "refGeneParser.hpp"
#include "utr3Finder.hpp"
#include "utr3FinderFuzzy.hpp"
#include "bppPredictFuzzy.hpp"

extern "C" {
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/utils.h>
}


namespace fs = boost::filesystem;
namespace qi = boost::spirit::qi;


/**
 * The truth values obtained from this function are inversed since it is the same used by Utr3FinderFuzzy.
 * Using the vector index as position reference for the motif (hence we got 6 items in this vector).
 */
std::vector<BppPredictFuzzy::tvFunction> BppPredictFuzzy::motifPositionToTruthValue  = {
	{BppPredictFuzzy::tvFunction(0.0, 0.25, 2.00)},
	{BppPredictFuzzy::tvFunction(0.0, 0.25, 2.00)},
	{BppPredictFuzzy::tvFunction(0.0, 0.25, 2.00)},
	{BppPredictFuzzy::tvFunction(0.0, 0.25, 2.00)},
	{BppPredictFuzzy::tvFunction(0.0, 0.25, 2.00)},
	{BppPredictFuzzy::tvFunction(0.0, 0.25, 2.00)}
};


/**
 * Initializing the lookup for transcript properties.
 */
RefGeneParser BppPredictFuzzy::refGen = RefGeneParser(fs::path("ucsc_data/refGene.txt"));


/**
 * Initializing static FASTA index. For later: Should get a proper initialization by argument from constructor.
 */
seqan::FaiIndex BppPredictFuzzy::faiIndex = []() {
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	std::cout << "Initializing FASTA index...";
	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		std::cerr << "could not open index file for " << referenceGenome << std::endl;
	}
	std::cout << "DONE" << std::endl;
	return faiIndex;
}();


/**
 * Initializing static map utrBppMap from a file to avoid unnecessary folding.
 */
BppPredictFuzzy::bppVectorPerTranscriptMap BppPredictFuzzy::utrBppMap = []() {
	std::cout << "Initializing BPP map...";
	fs::path bppFile = "../perf_testing/utrBppPerTranscript.txt";
	std::ifstream inFile(bppFile.string());
	std::string line;
	std::string transcriptId;
	BppPredictFuzzy::utrBppVector tempVec;
	BppPredictFuzzy::bppVectorPerTranscriptMap tempMap;
	size_t lineCount = 0;
	while (std::getline(inFile, line)) {
		switch (lineCount) {
		case 0:
		{	

			qi::parse(line.begin(), line.end(), 
				">" >> +qi::char_, 
				transcriptId);
			//std::cerr << transcriptId << std::endl;
			lineCount++;
			break;
		}
		case 1: 
		{	
			qi::parse(line.begin(), line.end(),
				 +(qi::double_ >> " "),
				 tempVec);
			if (tempMap.find(transcriptId) == tempMap.end()) {
				tempMap[transcriptId] = tempVec;
			}
			transcriptId.clear();
			tempVec.clear();
			lineCount = 0;
			break;
		}
		}
	}
	inFile.close();
	std::cout << "DONE" << std::endl;
	/*
	for (auto & pair : tempMap) {
		std::cerr << pair.first << ": ";
		for (auto & item : pair.second) {
			std::cerr << item << ", ";
		}
		std::cerr << std::endl;
	}
	*/
	return tempMap;	
		
}();


/*
 * Constructor.
 */
BppPredictFuzzy::BppPredictFuzzy(const std::string & id): 
	txId(id)
{
	if (utrBppMap.find(this->txId) != utrBppMap.end()) {
		this->maxBpp = &utrBppMap.find(this->txId)->second;
	}
	this->startPrediction();	
}


/**
 * Constructor. Used for specificity calculation. Allowing randomized input.
 */
BppPredictFuzzy::BppPredictFuzzy(const std::string & id, 
	const std::string & utr,
	const std::string & offSeq):
	utrSeq(utr),
	offsetSeq(offSeq),
	txId(id)
{
	//TODO: Provide randomized utr and offsetSeq and have it folded
	this->startPrediction();
}

/*
 * Destructor.
 */
BppPredictFuzzy::~BppPredictFuzzy() {}


/*
 * Entry point.
 */
void BppPredictFuzzy::startPrediction() {	
	this->evaluatePotentialPas();
	if (this->maxBpp == nullptr) {
		this->foldUtr();
	}
	/*
	for (auto & res : this->utr3FinderRes) {
		std::cerr << res.truthValue << " ";
	}
	std::cerr << std::endl;
	*/
	this->calcBppTruthValue();
	/*
	for (auto & res : this->utr3FinderRes) {
		std::cerr << res.truthValue << " ";
	}
	std::cerr << std::endl;
	*/
}


/**
 * Folds the UTR sequence provided by the user.
 */
void BppPredictFuzzy::foldUtr() {
	utrBppVector tempBppVec = std::vector<double>(this->utrSeq.size(), 0.0);
	const char * cStringUtr = this->utrSeq.c_str();
	char * mfe_structure = static_cast<char*>(vrna_alloc(sizeof(char) * (strlen(cStringUtr) + 1)));
	vrna_fold_compound_t * foldStruct = vrna_fold_compound (cStringUtr, NULL, VRNA_OPTION_PF | VRNA_OPTION_MFE);
	double mfe = (double)vrna_mfe(foldStruct, mfe_structure);
	vrna_exp_params_rescale(foldStruct, &mfe);
	vrna_mfe(foldStruct, mfe_structure);
	//std::cerr << test << std::endl << mfe_structure << std::endl;
	vrna_pf(foldStruct, NULL);
	//std::cerr << "exp_matrices->length: " <<  test_fold->exp_matrices->length << std::endl;
	vrna_plist_t * pl = vrna_plist_from_probs(foldStruct, 1e-5);
	//std::cerr << "maxBpp.size: " << this->maxBpp.size() << std::endl;
	for (size_t k = 0; pl[k].i>0 || pl[k].j>0; k++) {
		//std::cerr << pl[k].i << ", " << pl[k].j << ": " << pl[k].p << std::endl;
		if (pl[k].p > tempBppVec[pl[k].i] || pl[k].p > tempBppVec[pl[k].j]) {
			tempBppVec[pl[k].i - 1] = pl[k].p;
			tempBppVec[pl[k].j - 1] = pl[k].p;
		}	
	}
	free(pl);
	free(mfe_structure);
	vrna_fold_compound_free(foldStruct);
	
	utrBppMap[this->txId] = tempBppVec;
	this->maxBpp = &utrBppMap.find(this->txId)->second;
	/*
		BOOST_FOREACH(double & d, this->maxBpp) {
			std::cerr << d << ", ";
		}
		std::cerr << std::endl;
	*/
}


/**
 * Evaluate potential PAS in the utr sequence.  
 */
void BppPredictFuzzy::evaluatePotentialPas() {
	//Utr3FinderFuzzy should only search forward strand to avoid unnecessary calculations
	const bool searchBackward = false;
	//Give me a map with all the transcripts and their properties
	const RefGeneProperties & txProps = this->refGen.getValueByKey(this->txId);
        //Was the user-specified transcript ID found?
	if (txProps == this->refGen.emptyRefGeneProperties) {
		std::cerr << "Could not find " << this->txId << std::endl;
		return;
	}
	const std::string & strand = txProps.strand;
	size_t utrStart = txProps.cdsEnd;
	size_t utrEnd = txProps.txEnd;
	if (strand == "-") {
		utrStart = txProps.txStart;
		utrEnd = txProps.cdsStart;
	}
	size_t idx = polar::utility::getFastaIndex(txProps.chr);
	//Retrieving the sequence after the transcript end and concatenating it with the UTR sequence
	//So Utr3FinderFuzzy is able to analyze the PAS close to the transcript end.
	seqan::CharString seq;
	if (strand == "+") {
		seqan::readRegion(seq, BppPredictFuzzy::faiIndex, idx, utrEnd, utrEnd + this->txOffset + 1);
	} else {
		seqan::readRegion(seq, BppPredictFuzzy::faiIndex, idx, utrStart - (this->txOffset + 1), utrStart);
		seqan::reverseComplement(seq);
	}
	seqan::toLower(seq);
	this->utrSeq = polar::utility::getUtrSequence(txProps, this->faiIndex);
	this->offsetSeq = std::string(seqan::toCString(seq));
	//std::cerr << "txEnd: " << txEnd << ", chr: " << txProps.chr << ", txStart: " << txStart << std::endl;
	//std::cerr << utrSeq << sequence << std::endl;	
	SeqStruct sStruct = {
		utrSeq + offsetSeq,
		boost::none,
		boost::none,
		boost::none,
		boost::none,
		boost::none,
		boost::none,
		boost::none
	};
	Utr3FinderFuzzy candidates(sStruct, searchBackward);
	this->utr3FinderRes = candidates.getPolyaMotifPos();
}


/**
 * Calculating truth value.
 */
void BppPredictFuzzy::calcBppTruthValue() {
	//TODO: evaluate BPPs and calculate TV, then incorporate them with the TV from Utr3FinderFuzzy
	BOOST_FOREACH(resultStruct & candidate, this->utr3FinderRes) {
		size_t & pos = candidate.pos;
		if (pos + 6 > this->maxBpp->size()) continue;
		double truthValue = 0.0;
		double sumTv = 0.0;
		//std::cerr << "maxBpp.size(): " << this->maxBpp->size();
		//std::cerr << " position: " << pos << ", maxBpp[pos]: " << this->maxBpp->operator[](pos) << std::endl;
		
		for (size_t i = 0; i < 6; i++) {
			truthValue = this->getTruthValue((*(this->maxBpp))[pos + i], i);
			sumTv += truthValue;
		}
		//std::cerr << "sumTv: " << sumTv << " sumTv/6: " << sumTv / 6 << std::endl;
		candidate.truthValue += sumTv / 6;
	}
}


/**
 * Return truth value for a certain base pair probability and position.
 */
double BppPredictFuzzy::getTruthValue(const double & bpp, const size_t pos) const {
	const Utr3FinderFuzzy::UracilContent & bppTv = this->BppPredictFuzzy::motifPositionToTruthValue[pos];
	Utr3FinderFuzzy::UracilContent::straight interStraight = bppTv.getStraight();
	double uB = bppTv.getUpperBound();
	double lB = bppTv.getLowerBound();
	double maxTv = bppTv.getMaxTruthValue();
	
	if (bpp >= uB) {
		return 0.0;
	} else if (bpp <= lB) {
		return maxTv;
	} else {
		return maxTv - (interStraight.first * bpp + interStraight.second);
	}
}


/**
 * Return a vector of hits with a score (truth value) higher than a certain threshold.
 */
std::vector<BppPredictFuzzy::resultStruct> BppPredictFuzzy::getResults(const double & threshold) {
	std::vector<BppPredictFuzzy::resultStruct> resVec;
	BOOST_FOREACH(BppPredictFuzzy::resultStruct & candidate, this->utr3FinderRes) {
		if (candidate.truthValue >= threshold) {
			resVec.push_back(candidate);
		}
	}
	return resVec;
}


/**
 * Getter.
 */
std::string BppPredictFuzzy::getUtrSeq() const {
	return this->utrSeq;
}


/**
 * Getter.
 */
std::string BppPredictFuzzy::getOffsetSeq() const {
	return this->offsetSeq;
}


/**
 * Getter.
 */
std::vector<double> BppPredictFuzzy::getMaxBppVector() const {
	return *(this->maxBpp);
}


/**
 * Getter.
 */
std::string BppPredictFuzzy::getTxId() const {
	return this->txId;
}


#include <string>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
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


RefGeneParser BppPredictFuzzy::refGen = RefGeneParser(fs::path("ucsc_data/refGene.txt"));


/*
 * Constructor.
 */
BppPredictFuzzy::BppPredictFuzzy(const std::string & txId, const std::string & seq) {
	this->txId = txId;
	this->utrSeq = seq;
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
	this->foldUtr();
	this->evaluatePotentialPas();
}


/**
 * Folds the UTR sequence provided by the user.
 */
void BppPredictFuzzy::foldUtr() {
	const char * test = this->utrSeq.c_str();
	char * mfe_structure = static_cast<char*>(vrna_alloc(sizeof(char) * (strlen(test) + 1)));
	vrna_fold_compound_t * test_fold = vrna_fold_compound (test, NULL, VRNA_OPTION_PF | VRNA_OPTION_MFE);
	vrna_mfe(test_fold, mfe_structure);
	std::cerr << test << std::endl << mfe_structure << std::endl;
	vrna_pf(test_fold, NULL);
	std::cerr << "exp_matrices->length: " <<  test_fold->exp_matrices->length << std::endl;
	vrna_plist_t * pl = vrna_plist_from_probs(test_fold, 1e-5);
	
	for (size_t k = 0; pl[k].i>0 || pl[k].j>0; k++) {
		std::cerr << pl[k].i << ", " << pl[k].j << ": " << pl[k].p << std::endl;
	}
	
}


/**
 * Evaluate potential PAS in the utr sequence.
 */
void BppPredictFuzzy::evaluatePotentialPas() {
	const bool searchBackward = false;
	const RefGeneProperties & txProps = this->refGen.getValueByKey(this->txId);	
	const size_t cdsEnd = txProps.cdsEnd;
	const size_t txEnd = txProps.exonEnds.back();
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
		
	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		std::cerr << "could not open index file for " << referenceGenome << std::endl;
	}
	
	size_t idx = polar::utility::getFastaIndex(txProps.chr);

	seqan::CharString seq;	
	seqan::readRegion(seq, faiIndex, idx, cdsEnd, txEnd + this->txOffset);
	std::string sequence(seqan::toCString(seq));

	SeqStruct sStruct = {
		sequence,
		boost::none,
		boost::none,
		boost::none,
		boost::none,
		boost::none,
		boost::none,
		boost::none
	};

	Utr3FinderFuzzy candidates(sStruct, searchBackward);

	std::vector<Utr3Finder::Utr3FinderResult> results = candidates.getPolyaMotifPos();	
}





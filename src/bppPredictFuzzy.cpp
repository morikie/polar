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
	const char * cStringUtr = this->utrSeq.c_str();
	char * mfe_structure = static_cast<char*>(vrna_alloc(sizeof(char) * (strlen(cStringUtr) + 1)));
	vrna_fold_compound_t * foldStruct = vrna_fold_compound (cStringUtr, NULL, VRNA_OPTION_PF | VRNA_OPTION_MFE);
	vrna_mfe(foldStruct, mfe_structure);
	//std::cerr << test << std::endl << mfe_structure << std::endl;
	vrna_pf(foldStruct, NULL);
	//std::cerr << "exp_matrices->length: " <<  test_fold->exp_matrices->length << std::endl;
	vrna_plist_t * pl = vrna_plist_from_probs(foldStruct, 1e-5);
	
	for (size_t k = 0; pl[k].i>0 || pl[k].j>0; k++) {
		std::cerr << pl[k].i << ", " << pl[k].j << ": " << pl[k].p << std::endl;
	}
	
}


/**
 * Evaluate potential PAS in the utr sequence.  
 */
void BppPredictFuzzy::evaluatePotentialPas() {
	//Utr3FinderFuzzy should only search forward strand to avoid unnecessary calculations
	const bool searchBackward = false;
	//Give me a map with all the transcripts and their properties
	const RefGeneProperties & txProps = this->refGen.getValueByKey(this->txId);
        //Was the user specified transcript ID found?
	if (txProps == this->refGen.emptyRefGeneProperties) {
		std::cerr << "Could not find " << this->txId << std::endl;
		return;
	}
	const size_t & cdsEnd = txProps.cdsEnd;
	const size_t & txEnd = txProps.exonEnds.back();
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
		
	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		std::cerr << "could not open index file for " << referenceGenome << std::endl;
	}
	size_t idx = polar::utility::getFastaIndex(txProps.chr);
	
	//Retrieving the sequence after the transcript end and concatenating it with the UTR sequence
	//TODO: Consideration for "-"-strand transcripts!
	seqan::CharString seq;	
	seqan::readRegion(seq, faiIndex, idx, txEnd, this->txOffset);
	std::string sequence(seqan::toCString(seq));
	sequence = this->utrSeq + sequence;

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
	this->utr3FinderRes = candidates.getPolyaMotifPos();
	
}





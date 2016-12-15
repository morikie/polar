#ifndef __BPPPREDICTFUZZY_HPP__
#define __BPPPREDICTFUZZY_HPP__

#include <unordered_map>
#include <boost/filesystem/path.hpp>
#include <seqan/seq_io.h>
#include "../src/refGeneParser.hpp"
#include "utr3Finder.hpp"
#include "utr3FinderFuzzy.hpp"

extern "C" {
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/utils.h>
}

namespace fs = boost::filesystem;


class BppPredictFuzzy {
public:
	typedef std::string transcriptId;
	typedef Utr3FinderFuzzy::UracilContent tvFunction;
	typedef Utr3Finder::Utr3FinderResult resultStruct;
	typedef std::vector<double> utrBppVector;
	typedef std::unordered_map<transcriptId, utrBppVector> bppVectorPerTranscriptMap;
	//evaluation function; index of vector is reference to motif position
	static std::vector<tvFunction> motifPositionToTruthValue;
	//map that contains already calculated BPPs (read from a file)
	//If there is no entry in the map, the program will attempt to fold the sequence (might take long)
	static bppVectorPerTranscriptMap utrBppMap;
	//Fasta index needed by seqan's readRegion function to extract DNA sequences from a reference genome
	static seqan::FaiIndex faiIndex;

private:
	std::string utrSeq;
	std::string offsetSeq;
	std::string txId;
	static RefGeneParser refGen;
	//length of the considered bases after transcript end
	const size_t txOffset = 200;
	std::vector<resultStruct> utr3FinderRes;
	utrBppVector * maxBpp = nullptr;
	
public:
	BppPredictFuzzy(const std::string & id);
	BppPredictFuzzy(const std::string & id, const std::string & utr, const std::string & offSeq);
	~BppPredictFuzzy();
	
	std::vector<resultStruct> getResults(const double & threshold);
	std::vector<double> getMaxBppVector() const;
	std::string getUtrSeq() const;
	std::string getOffsetSeq() const;
	std::string getTxId() const;

private:
	void startPrediction();
	void foldUtr();
	void evaluatePotentialPas();
	void calcBppTruthValue();
	double getTruthValue(const double & bpp, const size_t pos) const;
	std::vector<resultStruct> getResults(const double & threshold) const;
};	


#endif /* __BPPPREDICTFUZZY_HPP__ */


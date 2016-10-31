#ifndef __BPPPREDICTFUZZY_HPP__
#define __BPPPREDICTFUZZY_HPP__

#include <boost/filesystem/path.hpp>
#include "../src/refGeneParser.hpp"
#include "utr3FinderFuzzy.hpp"

extern "C" {
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/utils.h>
}


class BppPredictFuzzy {
private:
	std::string utrSeq;
	std::string txId;
	static RefGeneParser refGen;
	const size_t txOffset = 200;

public:
	BppPredictFuzzy(const std::string & txId, const std::string & seq);
	~BppPredictFuzzy();

private:
	void startPrediction();
	void foldUtr();
	void evaluatePotentialPas();
};	


#endif /* __BPPPREDICTFUZZY_HPP__ */

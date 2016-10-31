#include <iostream>
#include <string>
#include "../src/bppPredictFuzzy.hpp"
#include "../src/refGeneParser.hpp"
#include "perf_bppFuzzy.hpp"

int main (int argc, char * argv[]) {
	std::string test = "gattaca";
	std::string txTestId = "NM_000123151";
	BppPredictFuzzy bppFuzzyTest = BppPredictFuzzy(txTestId, test);
}


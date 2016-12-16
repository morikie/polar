#include <algorithm>
#include <iostream>
#include <list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <boost/foreach.hpp>
#include "polarUtility.hpp"
#include "seqStruct.hpp"
#include "utr3FinderFuzzy.hpp"


///**
// * Initialization of static member that holds values for the ranges of the DSE location.
// * Taken from the original paper ("Prediction of non-canonical polyadenylation signals..."; doi: 10.1016/j.jbiosc.2009.01.001)
// */
//Utr3FinderFuzzy::pasToDseLocMap Utr3FinderFuzzy::dseLocMap = { 
//	{std::string("aataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 25), std::make_pair(35, 55))},
//	{std::string("attaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 25), std::make_pair(33, 60))},
//	{std::string("tataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 30), std::make_pair(40, 60))},
//	{std::string("agtaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 27), std::make_pair(36, 50))},
//	{std::string("aagaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 27), std::make_pair(35, 70))},
//	{std::string("aatata"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 22), std::make_pair(36, 70))},
//	{std::string("aataca"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
//	{std::string("cataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 27), std::make_pair(37, 55))},
//	{std::string("gataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
//	{std::string("aatgaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 25), std::make_pair(40, 65))},
//	{std::string("actaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 27), std::make_pair(37, 50))},
//	{std::string("aataga"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 30), std::make_pair(42, 60))}
//};


//Utr3FinderFuzzy::pasToUcontentMap Utr3FinderFuzzy::dseUracilMap = {
//	{std::string("aataaa"), Utr3FinderFuzzy::UracilContent(0.30, 0.80, 1.0)},
//	{std::string("attaaa"), Utr3FinderFuzzy::UracilContent(0.30, 0.80, 1.0)},
//	{std::string("tataaa"), Utr3FinderFuzzy::UracilContent(0.30, 0.80, 1.0)},
//	{std::string("agtaaa"), Utr3FinderFuzzy::UracilContent(0.30, 0.80, 1.0)},
//	{std::string("aagaaa"), Utr3FinderFuzzy::UracilContent(0.30, 0.80, 1.0)},
//	{std::string("aatata"), Utr3FinderFuzzy::UracilContent(0.30, 0.80, 1.0)},
//	{std::string("aataca"), Utr3FinderFuzzy::UracilContent(0.30, 0.80, 1.0)},
//	{std::string("cataaa"), Utr3FinderFuzzy::UracilContent(0.30, 0.80, 1.0)},
//	{std::string("gataaa"), Utr3FinderFuzzy::UracilContent(0.30, 0.80, 1.0)},
//	{std::string("aatgaa"), Utr3FinderFuzzy::UracilContent(0.30, 0.80, 1.0)},
//	{std::string("actaaa"), Utr3FinderFuzzy::UracilContent(0.30, 0.80, 1.0)},
//	{std::string("aataga"), Utr3FinderFuzzy::UracilContent(0.30, 0.80, 1.0)}
//};


/**
 * Initialization of static member that holds values for the lower/upper bound of the USE uracil content.
 * Taken from the original paper ("Prediction of non-canonical polyadenylation signals..."; doi: 10.1016/j.jbiosc.2009.01.001)
 */
//Utr3FinderFuzzy::pasToUcontentMap Utr3FinderFuzzy::useUracilMap = {
//	{std::string("aataaa"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.4)},
//	{std::string("attaaa"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.4)},
//	{std::string("tataaa"), Utr3FinderFuzzy::UracilContent(0.56, 0.66, 0.5)},
//	{std::string("agtaaa"), Utr3FinderFuzzy::UracilContent(0.56, 0.66, 0.2)},
//	{std::string("aagaaa"), Utr3FinderFuzzy::UracilContent(0.56, 0.66, 0.4)},
//	{std::string("aatata"), Utr3FinderFuzzy::UracilContent(0.56, 0.78, 0.4)},
//	{std::string("aataca"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.4)},
//	{std::string("cataaa"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.2)},
//	{std::string("gataaa"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.2)},
//	{std::string("aatgaa"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.2)},
//	{std::string("actaaa"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.4)},
//	{std::string("aataga"), Utr3FinderFuzzy::UracilContent(0.56, 0.66, 0.4)}
//};


/**
 * Initialization of static member that holds values for the threshold to determine an authentic PAS.
 * Taken from the original paper ("Prediction of non-canonical polyadenylation signals..."; doi: 10.1016/j.jbiosc.2009.01.001)
 */
//std::unordered_map<Utr3FinderFuzzy::motifSequence, double> Utr3FinderFuzzy::thresholdMap = {
//	{std::string("aataaa"), 0.51},
//	{std::string("attaaa"), 0.51},
//	{std::string("tataaa"), 0.74},
//	{std::string("agtaaa"), 0.51},
//	{std::string("aagaaa"), 0.51},
//	{std::string("aatata"), 0.74},
//	{std::string("aataca"), 0.55},
//	{std::string("cataaa"), 0.51},
//	{std::string("gataaa"), 0.51},
//	{std::string("aatgaa"), 0.51},
//	{std::string("actaaa"), 0.51},
//	{std::string("aataga"), 0.74}
//};


/**
 * Different TV for the dse location map.
 */
Utr3FinderFuzzy::pasToDseLocMap Utr3FinderFuzzy::dseLocMap = { 
	{std::string("aataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(13, 25), std::make_pair(35, 50))},
	{std::string("attaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(13, 25), std::make_pair(35, 50))},
	{std::string("tataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(13, 25), std::make_pair(35, 50))},
	{std::string("agtaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(13, 25), std::make_pair(35, 50))},
	{std::string("aagaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(13, 25), std::make_pair(35, 50))},
	{std::string("aatata"), Utr3FinderFuzzy::DseLocation(std::make_pair(13, 25), std::make_pair(35, 50))},
	{std::string("aataca"), Utr3FinderFuzzy::DseLocation(std::make_pair(13, 25), std::make_pair(35, 50))},
	{std::string("cataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(13, 25), std::make_pair(35, 50))},
	{std::string("gataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(13, 25), std::make_pair(35, 50))},
	{std::string("aatgaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(13, 25), std::make_pair(35, 50))},
	{std::string("actaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(13, 25), std::make_pair(35, 50))},
	{std::string("aataga"), Utr3FinderFuzzy::DseLocation(std::make_pair(13, 25), std::make_pair(35, 50))}
};


/**
 * Different TV for the dse location map.
 */
Utr3FinderFuzzy::pasToDseLocMap Utr3FinderFuzzy::dseDinucleoLocMap = { 
	{std::string("aataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 15), std::make_pair(30, 50))},
	{std::string("attaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 15), std::make_pair(30, 50))},
	{std::string("tataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 15), std::make_pair(30, 50))},
	{std::string("agtaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 15), std::make_pair(30, 50))},
	{std::string("aagaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 15), std::make_pair(30, 50))},
	{std::string("aatata"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 15), std::make_pair(30, 50))},
	{std::string("aataca"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 15), std::make_pair(30, 50))},
	{std::string("cataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 15), std::make_pair(30, 50))},
	{std::string("gataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 15), std::make_pair(30, 50))},
	{std::string("aatgaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 15), std::make_pair(30, 50))},
	{std::string("actaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 15), std::make_pair(30, 50))},
	{std::string("aataga"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 15), std::make_pair(30, 50))}
};


Utr3FinderFuzzy::pasToUcontentMap Utr3FinderFuzzy::dseDinucleoMap = {
	{std::string("aataaa"), Utr3FinderFuzzy::UracilContent(0.20, 0.60, 1.0)},
	{std::string("attaaa"), Utr3FinderFuzzy::UracilContent(0.20, 0.60, 1.0)},
	{std::string("tataaa"), Utr3FinderFuzzy::UracilContent(0.20, 0.60, 1.0)},
	{std::string("agtaaa"), Utr3FinderFuzzy::UracilContent(0.20, 0.60, 1.0)},
	{std::string("aagaaa"), Utr3FinderFuzzy::UracilContent(0.20, 0.60, 1.0)},
	{std::string("aatata"), Utr3FinderFuzzy::UracilContent(0.20, 0.60, 1.0)},
	{std::string("aataca"), Utr3FinderFuzzy::UracilContent(0.20, 0.60, 1.0)},
	{std::string("cataaa"), Utr3FinderFuzzy::UracilContent(0.20, 0.60, 1.0)},
	{std::string("gataaa"), Utr3FinderFuzzy::UracilContent(0.20, 0.60, 1.0)},
	{std::string("aatgaa"), Utr3FinderFuzzy::UracilContent(0.20, 0.60, 1.0)},
	{std::string("actaaa"), Utr3FinderFuzzy::UracilContent(0.20, 0.60, 1.0)},
	{std::string("aataga"), Utr3FinderFuzzy::UracilContent(0.20, 0.60, 1.0)}
};


/**
 * Initialization of static member that holds value for the lower/upper bound of the DSE uracil content.
 * Taken from the original paper ("Prediction of non-canonical polyadenylation signals..."; doi: 10.1016/j.jbiosc.2009.01.001)
 */
Utr3FinderFuzzy::pasToUcontentMap Utr3FinderFuzzy::dseUracilMap = {
	{std::string("aataaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 1.0)},
	{std::string("attaaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 1.0)},
	{std::string("tataaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 1.0)},
	{std::string("agtaaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 1.0)},
	{std::string("aagaaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 1.0)},
	{std::string("aatata"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 1.0)},
	{std::string("aataca"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 1.0)},
	{std::string("cataaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 1.0)},
	{std::string("gataaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 1.0)},
	{std::string("aatgaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 1.0)},
	{std::string("actaaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 1.0)},
	{std::string("aataga"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 1.0)}
};


/**
 * Map with uracil content TV scoring for the short region after the PAS.
 */
Utr3FinderFuzzy::pasToUcontentMap Utr3FinderFuzzy::dseShortUracilMap = {
	{std::string("aataaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("attaaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("tataaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("agtaaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("aagaaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("aatata"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("aataca"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("cataaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("gataaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("aatgaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("actaaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("aataga"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)}
};


/**
 *
 */
Utr3FinderFuzzy::pasToUcontentMap Utr3FinderFuzzy::useUracilMap = {
	{std::string("aataaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("attaaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("tataaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("agtaaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("aagaaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("aatata"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("aataca"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("cataaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("gataaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("aatgaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("actaaa"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)},
	{std::string("aataga"), Utr3FinderFuzzy::UracilContent(0.40, 0.70, 0.3)}
};


/**
 *
 */
std::unordered_map<Utr3FinderFuzzy::motifSequence, double> Utr3FinderFuzzy::thresholdMap = {
	{std::string("aataaa"), 1.23},
	{std::string("attaaa"), 1.23},
	{std::string("tataaa"), 1.23},
	{std::string("agtaaa"), 1.23},
	{std::string("aagaaa"), 1.23},
	{std::string("aatata"), 1.23},
	{std::string("aataca"), 1.23},
	{std::string("cataaa"), 1.23},
	{std::string("gataaa"), 1.23},
	{std::string("aatgaa"), 1.23},
	{std::string("actaaa"), 1.23},
	{std::string("aataga"), 1.23}
};


/**
 * Constructor.
 */
Utr3FinderFuzzy::Utr3FinderFuzzy(const SeqStruct & sSt, const bool & seaBack):
	Utr3Finder(sSt),
	searchBackward(seaBack)	
{
	if (this->searchBackward) {
		const std::string & seq = this->seqStruct.seq;
		std::transform(seq.rbegin(), seq.rend(), std::back_inserter(this->reversedSeq), polar::utility::complement);
	}
	this->findPolyaMotif();

}


/**
 * Destructor.
 */
Utr3FinderFuzzy::~Utr3FinderFuzzy() {}


/**
 * Predicts the Poly(A) motives in a sequence.
 */
void Utr3FinderFuzzy::findPolyaMotif() {
	const std::string & seq = this->seqStruct.seq;
	const std::string & revSeq = this->reversedSeq;
	std::vector<size_t> candidatePositions;
	std::vector<size_t> revCandidatePositions;
	//searching forward strand
	BOOST_FOREACH (const Utr3FinderFuzzy::pasToDseLocMap::value_type & v, Utr3FinderFuzzy::dseLocMap) {
		auto posIt = seq.begin();
		while ((posIt = std::search(posIt, seq.end(), v.first.begin(), v.first.end())) != seq.end()) {
			candidatePositions.push_back(std::distance(seq.begin(), posIt));
			posIt++;
		}
	}
	
	//searching backward strand
	if (this->searchBackward) {
		BOOST_FOREACH (const Utr3FinderFuzzy::pasToDseLocMap::value_type & v, Utr3FinderFuzzy::dseLocMap) {
			auto posIt = revSeq.begin();
			while ((posIt = std::search(posIt, revSeq.end(), v.first.begin(), v.first.end())) != revSeq.end()) {
				revCandidatePositions.push_back(std::distance(revSeq.begin(), posIt));
				posIt++;
			}
		}
	}
	//verifying forward candidates
	BOOST_FOREACH(const size_t & candPos, candidatePositions) {
		std::string motif(seq.begin() + candPos, seq.begin() + candPos + 6);
		double combDseTValue;
		double dseShortTValue;
		double useTvalue;
		double combDseDinucleoTValue;

		combDseTValue = calcCombinedDseTvalue(candPos, seq);
		dseShortTValue = calcDseShortTvalue(candPos, seq);
		useTvalue = calcUseTvalue(candPos, seq);
		combDseDinucleoTValue = calcCombinedDinucleoDseTvalue(candPos, seq);
		double finalTvalue = combDseTValue + combDseDinucleoTValue + useTvalue + dseShortTValue;
		if (finalTvalue >= thresholdMap.find(motif)->second) {
			polyaPosVector.push_back(Utr3FinderResult{
				candPos, 
				finalTvalue,
				(this->searchBackward) ? "+" : "0"
				});
		}
	}
	//verifying backward candidates
	if (this->searchBackward) {
		BOOST_FOREACH(const size_t & candPos, revCandidatePositions) {
			std::string motif(revSeq.begin() + candPos, revSeq.begin() + candPos + 6);
			double combDseTValue;
			double dseShortTValue;
			double useTvalue;
			double combDseDinucleoTValue;

			combDseTValue = calcCombinedDseTvalue(candPos, revSeq);
			dseShortTValue = calcDseShortTvalue(candPos, revSeq);
			useTvalue = calcUseTvalue(candPos, revSeq);
			combDseDinucleoTValue = calcCombinedDinucleoDseTvalue(candPos, revSeq);
			double finalTvalue = combDseTValue + combDseDinucleoTValue + useTvalue + dseShortTValue;
			if (finalTvalue >= thresholdMap.find(motif)->second) {
				polyaPosVector.push_back(Utr3FinderResult{
					seq.size() - 1 - candPos,
					finalTvalue,
					"-"
					});
			}
		}
	}
}


/*
 * Scans downstream of a PAS candidate(i.e. its position)  for an uracil-rich region (the DSE).
 * Returns the combined truth value of uracil content and location for a certain window.
 * 
 */
double Utr3FinderFuzzy::calcCombinedDseTvalue(const size_t & pos, const std::string & seq) {		
	std::string motif(seq.begin() + pos, seq.begin() + pos + 6);
	DseLocation & dseLoc = Utr3FinderFuzzy::dseLocMap.find(motif)->second;
		
	auto start = seq.begin() + pos + dseLoc.getLeftRange().first + 6;
	std::list<std::string::value_type> slidingWindow(start, start + this->windowSize);
	
	double dseLocTvalue = 0.0;
	double dseUcontentTvalue = 0.0;
	double truthValue = 0.0;
	double maxTruthValue = 0.0;
	size_t uracilCounter = std::count(slidingWindow.begin(),slidingWindow.end(), 't');
	
	std::string::const_iterator end;
	
	if (std::distance(start, seq.end()) <static_cast<int>(this->windowSize + 1)) return 0.0;
	
	if (dseLoc.getRightRange().second - dseLoc.getLeftRange().first < static_cast<size_t>(std::distance(start, seq.end()))) {
		end = start + dseLoc.getRightRange().second - dseLoc.getLeftRange().first;
	} else {
		end = seq.end() - (this->windowSize);
	}

	for (auto posIt = start; posIt != end; posIt++) {
		double uContent = static_cast<double>(uracilCounter) / this->windowSize;
		size_t distanceToPas = std::distance(start - dseLoc.getLeftRange().first, posIt);
		dseLocTvalue = locationTvalue(motif, distanceToPas, this->dseLocMap);
		dseUcontentTvalue = nucleotideContentTvalue(motif, uContent, this->dseUracilMap);
		truthValue = std::min(dseLocTvalue, dseUcontentTvalue);
		if (truthValue > maxTruthValue) maxTruthValue = truthValue;
		
		if (slidingWindow.front() == 't' && uracilCounter > 0) uracilCounter--;
		slidingWindow.pop_front();
		slidingWindow.push_back(*(posIt + this->windowSize));
		if (slidingWindow.back() == 't') uracilCounter++;
		
	}
	return maxTruthValue;
}


/**
 * Scans upstream of a potential PAS for a uracil-rich region.
 * Returns the max truth value.
 */
double Utr3FinderFuzzy::calcDseShortTvalue(const size_t & pos, const std::string & seq) {
	size_t searchRange = 10;
	size_t windowSize = 6;
	std::string motif(seq.begin() + pos, seq.begin() + pos + 6);
	std::string::const_iterator start = seq.begin() + pos + 6;

	std::list<std::string::value_type> slidingWindow(start, start + windowSize);
	
	double dseShortUcontentTvalue = 0.0;
	double maxTvalue = 0.0;
	size_t uracilCounter = std::count(slidingWindow.begin(), slidingWindow.end(), 't');

	if (static_cast<size_t>(std::distance(start, seq.end())) < windowSize + 1) return 0.0;
	
	std::string::const_iterator end;
	if (static_cast<size_t>(std::distance(start, seq.end())) < searchRange) {
		end = seq.end() - (windowSize);
	} else {
		end = start + searchRange - windowSize;
	}

	for (auto posIt = start; posIt != end; posIt++) {
		double uContent = static_cast<double>(uracilCounter) / windowSize;
		dseShortUcontentTvalue = nucleotideContentTvalue(motif, uContent, this->dseShortUracilMap);
		if (dseShortUcontentTvalue > maxTvalue) maxTvalue = dseShortUcontentTvalue;
		
		if (slidingWindow.front() == 't' && uracilCounter > 0) uracilCounter--;
		slidingWindow.pop_front();
		slidingWindow.push_back(*(posIt + windowSize));
		if (slidingWindow.back() == 't') uracilCounter++;
	}
	return maxTvalue;
}


/**
 * Scans upstream of a potential PAS for a uracil-rich region.
 * Returns the max truth value.
 */
double Utr3FinderFuzzy::calcUseTvalue(const size_t & pos, const std::string & seq) {
	size_t searchRange = 20;
	std::string motif(seq.begin() + pos, seq.begin() + pos + 6);
	std::string::const_reverse_iterator start = seq.rend() - pos;
	std::list<std::string::value_type> slidingWindow(start, start + this->windowSize);
	
	double useUcontentTvalue = 0.0;
	double maxTvalue = 0.0;
	size_t uracilCounter = std::count(slidingWindow.begin(), slidingWindow.end(), 't');

	if (static_cast<size_t>(std::distance(start, seq.rend())) < this->windowSize + 1) return 0.0;
	
	std::string::const_reverse_iterator end;
	if (static_cast<size_t>(std::distance(start, seq.rend())) < searchRange) {
		end = seq.rend() - (this->windowSize);
	} else {
		end = start + searchRange - (this->windowSize);
	}

	for (auto posIt = start; posIt != end; posIt++) {
		double uContent = static_cast<double>(uracilCounter) / this->windowSize;
		useUcontentTvalue = nucleotideContentTvalue(motif, uContent, this->useUracilMap);
		if (useUcontentTvalue > maxTvalue) maxTvalue = useUcontentTvalue;
		
		if (slidingWindow.front() == 't' && uracilCounter > 0) uracilCounter--;
		slidingWindow.pop_front();
		slidingWindow.push_back(*(posIt + this->windowSize));
		if (slidingWindow.back() == 't') uracilCounter++;
	}
	return maxTvalue;
}


/**
 * Calculates and returns the maximum dinucleotide concentration in a certain window.
 */
double Utr3FinderFuzzy::calcCombinedDinucleoDseTvalue(const size_t & pos, const std::string & seq) {
	size_t windowSize = 10;
	std::string dinucleotide("tg");
	std::string motif(seq.begin() + pos, seq.begin() + pos + 6);
	DseLocation & dseLoc = Utr3FinderFuzzy::dseDinucleoLocMap.find(motif)->second;
	auto start = seq.begin() + pos + 6 + dseLoc.getLeftRange().first;
	
	double dseDinucleoLocTvalue = 0.0;
	double dseUgContentTvalue = 0.0;
	double truthValue = 0.0;
	double maxTruthValue = 0.0;
	size_t ugCounter = 0;
	auto found = start;
	while ((found = std::search(found, start + windowSize, dinucleotide.begin(), dinucleotide.end())) != start + windowSize) 
		++found , ++ugCounter;
	
	std::string::const_iterator end;
	if (std::distance(start, seq.end()) < static_cast<int>(windowSize + 1)) return 0.0;
	
	//Use the whole considerable search frame if there's enough sequence available
	if (dseLoc.getRightRange().second - dseLoc.getLeftRange().first < static_cast<size_t>(std::distance(start, seq.end()))) {
		end = start + dseLoc.getRightRange().second - dseLoc.getLeftRange().first - windowSize;
	//otherwise use the remaining sequence
	} else {
		end = seq.end() - (windowSize);
	}
	
	//std::cout << pos << "  " << std::distance(start, end) << "  " << std::distance(start, seq.end()) << std::endl;
	//std::cerr << std::string(start - 6, seq.end()) << std::endl;
	//std::cerr << pos << std::endl;
	for (auto posIt = start; posIt != end; ++posIt) {
		double ugContent = static_cast<double>(ugCounter) / (windowSize / 2);
		size_t distanceToPas = std::distance(start - dseLoc.getLeftRange().first, posIt);
		dseDinucleoLocTvalue = locationTvalue(motif, distanceToPas, this->dseDinucleoLocMap);
		dseUgContentTvalue = nucleotideContentTvalue(motif, ugContent, this->dseDinucleoMap);
		truthValue = std::min(dseDinucleoLocTvalue, dseUgContentTvalue);
		if (truthValue > maxTruthValue) maxTruthValue = truthValue;
		
		//check if dinucleotide leaves the window (decrement)
		if (seq.compare(static_cast<size_t>(std::distance(seq.begin(), posIt)), 2, dinucleotide) == 0 && ugCounter > 0)
			--ugCounter;
		
		//check if dinucleotide enters the window (increment)
		//std::cout << static_cast<size_t>(std::distance(seq.begin(), posIt)) + windowSize - 1 << std::endl;
		if (seq.compare(static_cast<size_t>(std::distance(seq.begin(), posIt)) + windowSize, 2, dinucleotide) == 0)
			++ugCounter;
		
		//std::cerr << std::string(posIt, posIt + 10) << ", ugCounter: " << ugCounter << std::endl;
		

	}
	return maxTruthValue;
}


/**
 * Set the threshold map to determine a correct prediction.
 */
void Utr3FinderFuzzy::setThresholdMap(std::unordered_map<std::string, double> & map) {
	Utr3FinderFuzzy::thresholdMap = map;
}


/**
 * Checks if a variant hits the Poly(A) motif.
 */
bool Utr3FinderFuzzy::isMutationInMotif() const {
	return false;
}


/**
 * Returns the sequence that was searched on.
 */
std::string Utr3FinderFuzzy::getSequence() const {
	return this->seqStruct.seq;
}


/**
 * Returns the reverse complemented sequence that was searched on.
 */
std::string Utr3FinderFuzzy::getRevComplementSeq() const {
	return this->reversedSeq;
}


/**
 * Returns all authentic PAS found in the sequence.
 * [Probably needs to be deleted. Redundant/not needed?]
 */
std::string Utr3FinderFuzzy::getMotifSequence(const Utr3FinderResult & result) const {
	if (result.pos == Utr3Finder::noHitPos || result.pos > this->seqStruct.seq.size() - 6) return std::string();
	if (result.strand =="+") {
		auto motifStart = this->seqStruct.seq.begin() + result.pos;
		auto motifEnd = this->seqStruct.seq.begin() + result.pos + 6;
		return std::string(motifStart, motifEnd);	
	} else if (result.strand == "-") {
		auto motifStart = this->reversedSeq.begin() + result.pos;
		auto motifEnd = this->reversedSeq.begin() + result.pos + 6;
		return std::string(motifStart, motifEnd);	
	} else {
		return std::string();
	}
}


/**
 * Returns the positions of the found PAS.
 */
std::vector<Utr3Finder::Utr3FinderResult> Utr3FinderFuzzy::getPolyaMotifPos() const {
	return this->polyaPosVector;
}


/**
 * Writes information about the variant (from VCF e.g.) and the predicted PASs.
 */
void Utr3FinderFuzzy::writeInfo() const {
	std::stringstream ss;
	ss << "Poly(A) pos: ";
	if (this->polyaPosVector.empty()) {
		ss << "none found";
		return;
	}
	for (auto it = this->polyaPosVector.begin(); it != polyaPosVector.end(); it++) {
		ss << it->pos << "(" << it->strand << ") ";
	}
	if (this->seqStruct.genomicPos) ss << ", gPos: " << *(this->seqStruct.genomicPos);
	if (this->seqStruct.chrom) ss << ", chrom: " << *(this->seqStruct.chrom);
	
	std::cerr << ss.str() << std::endl;
}


/**
 * Returns the truth value for a given position of a potential DSE.
 */
double Utr3FinderFuzzy::locationTvalue(const std::string & pas, const size_t & pos, const Utr3FinderFuzzy::pasToDseLocMap & map) const {
	const DseLocation & dseLoc = map.find(pas)->second;
	Utr3FinderFuzzy::DseLocation::range leftRange = dseLoc.getLeftRange();
	Utr3FinderFuzzy::DseLocation::range rightRange = dseLoc.getRightRange();
	Utr3FinderFuzzy::DseLocation::straight leftStraight = dseLoc.getLeftStraight();
	Utr3FinderFuzzy::DseLocation::straight rightStraight = dseLoc.getRightStraight();

	if (pos >= leftRange.second && pos <= rightRange.first) {
		return 1.0;
	} else if (pos < leftRange.first || pos > rightRange.second) {
		return 0.0;
	} else if (pos > leftRange.first && pos < leftRange.second) {
		return leftStraight.first * static_cast<double>(pos) + leftStraight.second;
	} else {
		return rightStraight.first * static_cast<double>(pos) + rightStraight.second;
	}
}


/**
 * Returns the truth value for a given uracil content for a DSE.
 */
double Utr3FinderFuzzy::nucleotideContentTvalue(const std::string & pas, const double & uContent, const Utr3FinderFuzzy::pasToUcontentMap & map) const {	
	const UracilContent & nucleotideContent = map.find(pas)->second;
	Utr3FinderFuzzy::UracilContent::straight interStraight = nucleotideContent.getStraight();
	double uB = nucleotideContent.getUpperBound();
	double lB = nucleotideContent.getLowerBound();
	double maxTv = nucleotideContent.getMaxTruthValue();
	
	if (uContent >= uB) {
		return maxTv;
	} else if (uContent <= lB) {
		return 0.0;
	} else {
		return interStraight.first * uContent + interStraight.second;
	}
}


/**
 * Constructor DseLocation.
 */
Utr3FinderFuzzy::DseLocation::DseLocation(Utr3FinderFuzzy::DseLocation::range p, Utr3FinderFuzzy::DseLocation::range n):
	positiveIntermediate(p),
	negativeIntermediate(n)
{
	this->calcStraights();	
}


/**
 * Destructor
 */
Utr3FinderFuzzy::DseLocation::~DseLocation() {}


/**
 * Calculates the slopes and intercepts of the two straights (the two sides of the trapezoid) for the values between zero and one.
 */
void Utr3FinderFuzzy::DseLocation::calcStraights() {
	if (positiveIntermediate.second < positiveIntermediate.first) 
		throw std::invalid_argument("DseLocation: invalid range arguments");
	if (negativeIntermediate.second < negativeIntermediate.first) 
		throw std::invalid_argument("DseLocation: invalid range arguments");
	if (positiveIntermediate.first > negativeIntermediate.first || positiveIntermediate.first > negativeIntermediate.second) 
		throw std::invalid_argument("DseLocation: first range has larger values than second range");

	double slopePositiveIntermediate = 1.0 / 
		(static_cast<double>(positiveIntermediate.second) - static_cast<double>(positiveIntermediate.first));
	double slopeNegativeIntermediate = - 1.0 / 
		(static_cast<double>(negativeIntermediate.second) - static_cast<double>(negativeIntermediate.first));
	
	double interceptPositiveIntermediate = 1.0 - slopePositiveIntermediate * static_cast<double>(positiveIntermediate.second);
	double interceptNegativeIntermediate = 1.0 - slopeNegativeIntermediate * static_cast<double>(negativeIntermediate.first);

	this->positiveStraight = std::make_pair(slopePositiveIntermediate, interceptPositiveIntermediate);
	this->negativeStraight = std::make_pair(slopeNegativeIntermediate, interceptNegativeIntermediate);
}


/**
 * Getter.
 */
Utr3FinderFuzzy::DseLocation::range Utr3FinderFuzzy::DseLocation::getLeftRange() const {
	return this->positiveIntermediate;
}


/**
 * Getter.
 */
Utr3FinderFuzzy::DseLocation::range Utr3FinderFuzzy::DseLocation::getRightRange() const {
	return this->negativeIntermediate;
}


/**
 * Getter.
 */
Utr3FinderFuzzy::DseLocation::straight Utr3FinderFuzzy::DseLocation::getLeftStraight() const {
	return this->positiveStraight;
}


/**
 * Getter.
 */
Utr3FinderFuzzy::DseLocation::straight Utr3FinderFuzzy::DseLocation::getRightStraight() const  {
	return this->negativeStraight;
}


/**
 * Constructor UracilContent.
 */
Utr3FinderFuzzy::UracilContent::UracilContent(double lB, double uB, double mTv):
	lowerBound(lB),
	upperBound(uB),
	maxTruthValue(mTv) 
{
	this->calcStraight();	
}


/**
 * Destructor
 */
Utr3FinderFuzzy::UracilContent::~UracilContent() {}


/**
 * Calculates the slope and interecept of the straight for the values between zero and one.
 */
void Utr3FinderFuzzy::UracilContent::calcStraight() {
	if (upperBound < lowerBound) 
		throw std::invalid_argument("UracilContent: lower bound bigger than upper bound");
	if (this->maxTruthValue < 0.0) 
		throw std::invalid_argument("UracilContent: maximum truth value is lower than zero");
	
	double slopeIntermediate = this->maxTruthValue / (static_cast<double>(upperBound) - static_cast<double>(lowerBound));
	double interceptIntermediate = this->maxTruthValue - slopeIntermediate * static_cast<double>(upperBound);

	this->intermediate = std::make_pair(slopeIntermediate, interceptIntermediate);
}	


/**
 * Getter.
 */
double Utr3FinderFuzzy::UracilContent::getLowerBound() const {
	return this->lowerBound;
}


/**
 * Getter.
 */
double Utr3FinderFuzzy::UracilContent::getUpperBound() const {
	return this->upperBound;
}


/**
 * Getter.
 */
double Utr3FinderFuzzy::UracilContent::getMaxTruthValue() const {
	return this->maxTruthValue;
}


/**
 * Getter.
 */
Utr3FinderFuzzy::UracilContent::straight Utr3FinderFuzzy::UracilContent::getStraight() const {
	return this->intermediate;
}


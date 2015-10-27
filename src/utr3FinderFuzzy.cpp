#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include "seqStruct.hpp"
#include "utr3FinderFuzzy.hpp"


/**
 * Initialization of static member that holds values for the ranges of the DSE location.
 * Taken from the original paper ("Prediction of non-canonical polyadenylation signals..."; doi: 10.1016/j.jbiosc.2009.01.001)
 */
std::unordered_map<std::string, Utr3FinderFuzzy::DseLocation> Utr3FinderFuzzy::dseLocMap = { 
	{std::string("aataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 25), std::make_pair(35, 55))},
	{std::string("attaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 25), std::make_pair(33, 60))},
	{std::string("tataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 30), std::make_pair(40, 60))},
	{std::string("agtaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 27), std::make_pair(36, 50))},
	{std::string("aagaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 27), std::make_pair(35, 70))},
	{std::string("aatata"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 22), std::make_pair(36, 70))},
	{std::string("aataca"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
	{std::string("cataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 27), std::make_pair(37, 55))},
	{std::string("gataaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
	{std::string("aatgaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 25), std::make_pair(40, 65))},
	{std::string("actaaa"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 27), std::make_pair(37, 50))},
	{std::string("aataga"), Utr3FinderFuzzy::DseLocation(std::make_pair(10, 30), std::make_pair(42, 60))}
};


/**
 * Initialization of static member that holds value for the lower/upper bound of the DSE uracil content.
 * Taken from the original paper ("Prediction of non-canonical polyadenylation signals..."; doi: 10.1016/j.jbiosc.2009.01.001)
 */
std::unordered_map<std::string, Utr3FinderFuzzy::UracilContent> Utr3FinderFuzzy::dseUracilMap = {
	{std::string("aataaa"), Utr3FinderFuzzy::UracilContent(0.33, 0.78, 1.0)},
	{std::string("attaaa"), Utr3FinderFuzzy::UracilContent(0.33, 0.78, 1.0)},
	{std::string("tataaa"), Utr3FinderFuzzy::UracilContent(0.33, 0.78, 1.0)},
	{std::string("agtaaa"), Utr3FinderFuzzy::UracilContent(0.33, 0.78, 1.0)},
	{std::string("aagaaa"), Utr3FinderFuzzy::UracilContent(0.33, 0.78, 1.0)},
	{std::string("aatata"), Utr3FinderFuzzy::UracilContent(0.33, 0.78, 1.0)},
	{std::string("aataca"), Utr3FinderFuzzy::UracilContent(0.33, 0.78, 1.0)},
	{std::string("cataaa"), Utr3FinderFuzzy::UracilContent(0.33, 0.78, 1.0)},
	{std::string("gataaa"), Utr3FinderFuzzy::UracilContent(0.33, 0.78, 1.0)},
	{std::string("aatgaa"), Utr3FinderFuzzy::UracilContent(0.33, 0.78, 1.0)},
	{std::string("actaaa"), Utr3FinderFuzzy::UracilContent(0.33, 0.78, 1.0)},
	{std::string("aataga"), Utr3FinderFuzzy::UracilContent(0.33, 0.78, 1.0)}
};


/**
 * Initialization of static member that holds values for the lower/upper bound of the USE uracil content.
 * Taken from the original paper ("Prediction of non-canonical polyadenylation signals..."; doi: 10.1016/j.jbiosc.2009.01.001)
 */
std::unordered_map<std::string, Utr3FinderFuzzy::UracilContent> Utr3FinderFuzzy::useUracilMap = {
	{std::string("aataaa"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.4)},
	{std::string("attaaa"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.4)},
	{std::string("tataaa"), Utr3FinderFuzzy::UracilContent(0.56, 0.66, 0.5)},
	{std::string("agtaaa"), Utr3FinderFuzzy::UracilContent(0.56, 0.66, 0.2)},
	{std::string("aagaaa"), Utr3FinderFuzzy::UracilContent(0.56, 0.66, 0.4)},
	{std::string("aatata"), Utr3FinderFuzzy::UracilContent(0.56, 0.78, 0.4)},
	{std::string("aataca"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.4)},
	{std::string("cataaa"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.2)},
	{std::string("gataaa"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.2)},
	{std::string("aatgaa"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.2)},
	{std::string("actaaa"), Utr3FinderFuzzy::UracilContent(0.66, 0.78, 0.4)},
	{std::string("aataga"), Utr3FinderFuzzy::UracilContent(0.56, 0.66, 0.4)}
};


/**
 * Constructor.
 */
Utr3FinderFuzzy::Utr3FinderFuzzy(const SeqStruct & sSt):
	Utr3Finder(sSt)
{	
	this->findPolyaMotif();
}


/**
 * Destructor.
 */
Utr3FinderFuzzy::~Utr3FinderFuzzy() {}


/**
 * Predicts the Poly(A) motif in a sequence.
 */
void Utr3FinderFuzzy::findPolyaMotif() {
	std::string test("hello");	
	
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
	return std::string();
}


/**
 * Returns all authentic PAS found in the sequence.
 */
std::vector<std::string> Utr3FinderFuzzy::getMotifSequence() const {
	return std::vector<std::string>(1, std::string());
}


/**
 * Returns the positions of the found PAS.
 */
std::vector<size_t> Utr3FinderFuzzy::getPolyaMotifPos() const {
	return std::vector<size_t>(1, 0);
}


/**
 * Writes information about the variant (from VCF e.g.) and the predicted PASs.
 */
void Utr3FinderFuzzy::writeInfo() const {}


/**
 * Returns the truth value for a given position of a potential DSE.
 */
double Utr3FinderFuzzy::getDseLocationTvalue(std::string pas, size_t pos) const {
	Utr3FinderFuzzy::DseLocation::range leftRange = Utr3FinderFuzzy::dseLocMap.find(pas)->second.getLeftRange();
	Utr3FinderFuzzy::DseLocation::range rightRange = Utr3FinderFuzzy::dseLocMap.find(pas)->second.getRightRange();
	Utr3FinderFuzzy::DseLocation::straight leftStraight = Utr3FinderFuzzy::dseLocMap.find(pas)->second.getLeftStraight();
	Utr3FinderFuzzy::DseLocation::straight rightStraight = Utr3FinderFuzzy::dseLocMap.find(pas)->second.getRightStraight();

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
 * Returns the truth value for a given uracil content for a DSE (9nt window).
 */
double Utr3FinderFuzzy::getDseUcontentTvalue(std::string pas, double uContent) const {
	Utr3FinderFuzzy::UracilContent::straight interStraight = Utr3FinderFuzzy::dseUracilMap.find(pas)->second.getStraight();
	double uB = Utr3FinderFuzzy::dseUracilMap.find(pas)->second.getUpperBound();
	double lB = Utr3FinderFuzzy::dseUracilMap.find(pas)->second.getLowerBound();
	double maxTv = Utr3FinderFuzzy::dseUracilMap.find(pas)->second.getMaxTruthValue();
	
	if (uContent >= uB) {
		return maxTv;
	} else if (uContent <= lB) {
		return 0.0;
	} else {
		return interStraight.first * uContent + interStraight.second;
	}
}


/**
 * Returns the truth value for a given uracil content for a USE.
 */
double Utr3FinderFuzzy::getUseUcontentTvalue(std::string pas, double uContent) const {
	Utr3FinderFuzzy::UracilContent::straight interStraight = Utr3FinderFuzzy::useUracilMap.find(pas)->second.getStraight();
	double uB = Utr3FinderFuzzy::useUracilMap.find(pas)->second.getUpperBound();
	double lB = Utr3FinderFuzzy::useUracilMap.find(pas)->second.getLowerBound();
	double maxTv = Utr3FinderFuzzy::useUracilMap.find(pas)->second.getMaxTruthValue();
	
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
	double interceptNegativeIntermediate = 1.0 - slopeNegativeIntermediate * static_cast<double>(negativeIntermediate.second);

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


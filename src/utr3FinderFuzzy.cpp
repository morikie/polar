#include "seqStruct.hpp"
#include "utr3FinderFuzzy.hpp"


Utr3FinderFuzzy::Utr3FinderFuzzy(const SeqStruct & sSt):
	Utr3Finder(sSt)
{
	this->findPolyaMotif();
}


Utr3FinderFuzzy::~Utr3FinderFuzzy() {}

void Utr3FinderFuzzy::findPolyaMotif() {
/* implementation */
}


bool Utr3FinderFuzzy::isMutationInMotif() const {
	return false;
}


std::string Utr3FinderFuzzy::getSequence() const {
	return std::string();
}


std::string Utr3FinderFuzzy::getMotifSequence() const {
	return std::string();
}


size_t Utr3FinderFuzzy::getPolyaMotifPos() const {
	return 0;
}


void Utr3FinderFuzzy::writeInfo() const {}


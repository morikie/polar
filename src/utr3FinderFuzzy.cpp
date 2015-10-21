#include <string>
#include <vector>
#include "seqStruct.hpp"
#include "utr3FinderFuzzy.hpp"


class DseLocation {
	
};


class DseUracil {

};


class UseUracil {

};


Utr3FinderFuzzy::Utr3FinderFuzzy(const SeqStruct & sSt):
	Utr3Finder(sSt)
{
	this->findPolyaMotif();
}


Utr3FinderFuzzy::~Utr3FinderFuzzy() {}


void Utr3FinderFuzzy::findPolyaMotif() {
	
}


bool Utr3FinderFuzzy::isMutationInMotif() const {
	return false;
}


std::string Utr3FinderFuzzy::getSequence() const {
	return std::string();
}


std::vector<std::string> Utr3FinderFuzzy::getMotifSequence() const {
	return std::vector<std::string>(1, std::string());
}


std::vector<size_t> Utr3FinderFuzzy::getPolyaMotifPos() const {
	return std::vector<size_t>(1, 0);
}


void Utr3FinderFuzzy::writeInfo() const {}


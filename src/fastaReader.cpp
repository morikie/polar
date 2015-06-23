//fastaReader.cpp

#include <iostream>
#include "fastaReader.hpp"


namespace fs = boost::filesystem;


FastaReader::FastaReader(const fs::path & f) {
	this->file = f;	
	this->parse();
}


FastaReader::~FastaReader() {}


const fs::path & FastaReader::getFile() const {
	return this->file;
}


const FastaReader::fastaVector::const_iterator FastaReader::getBeginIterator() const {
	return this->fV.cbegin();
}


const FastaReader::fastaVector::const_iterator FastaReader::getEndIterator() const {
	return this->fV.cend();
}


void FastaReader::parse() {
	
}


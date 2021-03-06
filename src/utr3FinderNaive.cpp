#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include "utr3FinderNaive.hpp"


/**
 * Constructor.
 */
Utr3FinderNaive::Utr3FinderNaive(const SeqStruct & sSt):
	Utr3Finder(sSt)
{
	this->findPolyaMotif();
}


/**
 * Destructor.
 */
Utr3FinderNaive::~Utr3FinderNaive() {}


/**
 * Predicts the location of the UTR3's Poly(A) motif.
 */
void Utr3FinderNaive::findPolyaMotif() {
	//transcript length is not always the same as sequence size
	size_t seqLength;
	if (this->seqStruct.txLength) {
		seqLength = *(this->seqStruct.txLength);
	} else {
		seqLength = this->seqStruct.seq.size();
	}
	size_t hitPos = Utr3Finder::noHitPos;
	//length of the sequence searched from transcript end
	size_t searchRange = 100;
	//excepting remnants of the Poly(A) tail from the search
	size_t searchOffset = 0;
	auto it = this->seqStruct.seq.rbegin(); 
	while (*it == 'a') {
		searchOffset++;
		it++;
	}
	
	BOOST_FOREACH (const std::string & pattern, Utr3Finder::rHexamers) {
		//std::search returns an iterator to the start of the match
		//(in this case to the "end" of the match since we use reverse_iterators)
		auto posIt = std::search(this->seqStruct.seq.rbegin() + searchOffset, 
			this->seqStruct.seq.rbegin() + searchOffset + searchRange, 
			pattern.begin(), pattern.end());
		
		if (posIt != this->seqStruct.seq.rbegin() + searchOffset + searchRange) {
			//Since positions are indexed starting with zero we have to substract 1 
			//from the seqLength and the std::distance (+5 instead of +6) value
			hitPos = (seqLength - 1) - (std::distance(this->seqStruct.seq.rbegin(), posIt + 5));
			if (this->seqStruct.strand) { 
				this->polyaPosVector.push_back(
					Utr3Finder::Utr3FinderResult{hitPos, 1, *(this->seqStruct.strand)}
				);
			} else {	
				this->polyaPosVector.push_back(
					Utr3Finder::Utr3FinderResult{hitPos, 1, ""}
				);
			}
			break;
		} 	
	}
	if (hitPos == Utr3Finder::noHitPos && this->seqStruct.strand) {
		this->polyaPosVector.push_back(Utr3FinderResult{hitPos, 0, *(this->seqStruct.strand)});
	} else if (hitPos == Utr3Finder::noHitPos) {
		this->polyaPosVector.push_back(Utr3FinderResult{hitPos, 0, ""});
	}
}


/**
 * Checks if variant is in the Poly(A) motif.
 */
bool Utr3FinderNaive::isMutationInMotif() const {
	if (this->polyaPosVector[0].pos == Utr3Finder::noHitPos) return false;
	if (! this->seqStruct.mutation) return false;
	if (! this->seqStruct.utr3Start) return false;

	int diff = this->seqStruct.mutation->getMutPosition() + static_cast<int>(*this->seqStruct.utr3Start) -
		static_cast<int>(this->polyaPosVector[0].pos);
	if (diff < 6 && diff >= 0) {
		return true;
	} else {
		return false;
	}
}


/**
 * Returns the Poly(A) motif position in the transcript.
 */
std::vector<Utr3Finder::Utr3FinderResult> Utr3FinderNaive::getPolyaMotifPos() const {
	return this->polyaPosVector;
}


/**
 * Returns the motif hexamer.
 */
std::string Utr3FinderNaive::getMotifSequence(const Utr3FinderResult & result) const {
	if (result.pos == Utr3Finder::noHitPos || result.pos > this->seqStruct.seq.size() - 6) return std::string();

	auto motifStart = this->seqStruct.seq.begin() + result.pos;
	auto motifEnd = this->seqStruct.seq.begin() + result.pos + 6;
	return std::string(motifStart, motifEnd);
}


/**
 * Return string of the searched sequence.
 */
std::string Utr3FinderNaive::getSequence() const {
	return this->seqStruct.seq;
}


/**
 * Writes a string to the error console with information about the transcript and mutation (if available).
 */
void Utr3FinderNaive::writeInfo() const {
	std::stringstream ss;
	ss << "Poly(A) Pos: " << this->polyaPosVector[0].pos;
	if (this->seqStruct.utr3Start) ss << ", utr3Start: " << *(this->seqStruct.utr3Start);
	if (this->seqStruct.txLength) ss << ", txLength: " << *(this->seqStruct.txLength);
	if (this->seqStruct.chrom) ss << ", chrom: " <<  *(this->seqStruct.chrom);
	if (this->seqStruct.genomicPos) ss << ", gPos: " <<  *(this->seqStruct.genomicPos);
	if (this->seqStruct.seqId) ss << ", seqID: " << *(this->seqStruct.seqId);
	if (this->seqStruct.strand) ss << ", strand: " << *(this->seqStruct.strand);
	if (this->seqStruct.mutation && this->seqStruct.utr3Start) 
		ss << ", MutPos: " << *(this->seqStruct.utr3Start) + this->seqStruct.mutation->getMutPosition();
	
	std::cerr << ss.str() << std::endl;
}


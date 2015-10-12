#include <algorithm>
#include <iostream>
#include <sstream>
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
	size_t seqLength = seqStruct.txLength;
	//cutoff is difference between true length (like str.size())
	//of the sequence and the designated length from external sources.
	size_t cutoff = seqStruct.seq.size() - seqLength;
	size_t utr3Length = seqLength - this->seqStruct.utr3Start;
	size_t hitPos = Utr3FinderNaive::noHitPos;
	//length of the sequence searched from transcript end
	size_t searchRange = 100;
	//excepting remnants of the Poly(A) tail from the search
	size_t searchOffset = 0;
	auto it = this->seqStruct.seq.rbegin() + cutoff; 
	while (*it == 'a') {
		searchOffset++;
		it++;
	}
	if (searchRange > utr3Length) searchRange = utr3Length; 
	
	BOOST_FOREACH (const std::string & pattern, this->rHexamers) {
		//std::search returns an iterator to the start of the match
		//(in this case to the "end" of the match since we use reverse_iterators)
		auto posIt = std::search(this->seqStruct.seq.rbegin() + cutoff + searchOffset, 
			this->seqStruct.seq.rbegin() + cutoff + searchOffset + searchRange, 
			pattern.begin(), pattern.end());
		
		if (posIt != this->seqStruct.seq.rbegin() + cutoff + searchOffset + searchRange) {
			//Since positions are indexed starting with zero we have to substract 1 
			//from the seqLength and the std::distance (+5 instead of +6) value
			hitPos = (seqLength - 1) - (std::distance(this->seqStruct.seq.rbegin() + cutoff, posIt + 5));
			this->polyaMotifPos = hitPos;
			break;
		} 	
	}
	if (hitPos == Utr3FinderNaive::noHitPos) this->polyaMotifPos = hitPos;
}


/**
 * Checks if variant is in the Poly(A) motif.
 */
bool Utr3FinderNaive::isMutationInMotif() const {
	if (this->polyaMotifPos == Utr3FinderNaive::noHitPos) return false;
	if (! this->seqStruct.mutation) return false;

	int diff = this->seqStruct.mutation->getMutPosition() + static_cast<int>(this->seqStruct.utr3Start) -
		static_cast<int>(this->polyaMotifPos);
	if (diff < 6 && diff >= 0) {
		return true;
	} else {
		return false;
	}
}


/**
 * Returns the Poly(A) motif position in the transcript.
 */
size_t Utr3FinderNaive::getPolyaMotifPos() const {
	return this->polyaMotifPos;
}


/**
 * Returns the motif hexamer.
 */
std::string Utr3FinderNaive::getMotifSequence() const {
	if (this->polyaMotifPos == Utr3FinderNaive::noHitPos) return std::string();

	auto motifStart = this->seqStruct.seq.begin() + this->polyaMotifPos;
	auto motifEnd = this->seqStruct.seq.begin() + this->polyaMotifPos + 6;
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
	ss << "Poly(A) Pos: " << this->polyaMotifPos 
		<< ", utr3Start: " << this->seqStruct.utr3Start
		<< ", txLength: " << this->seqStruct.txLength;
	if (this->seqStruct.chrom) ss << ", chrom: " <<  *this->seqStruct.chrom;
	if (this->seqStruct.genomicPos) ss << ", gPos: " <<  *this->seqStruct.genomicPos;
	if (this->seqStruct.seqId) ss << ", seqID: " << *this->seqStruct.seqId;
	if (this->seqStruct.strand) ss << ", strand: " << *this->seqStruct.strand;
	if (this->seqStruct.mutation) ss << ", MutPos: " << this->seqStruct.utr3Start + this->seqStruct.mutation->getMutPosition();
	
	std::cerr << ss.str() << std::endl;
}



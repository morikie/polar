#include <algorithm>
#include <sstream>
#include <vector>
#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <seqan/find.h>
#include "utr3Finder.hpp"


/**
 * Vector of known Poly(A) cleavage motifs to be not pathogenic. Sorted by frequency (highest first).
 */
const std::vector<std::string> Utr3Finder::hexamers = {
	"aataaa",
	"attaaa",
	"tataaa",
	"agtaaa",
	"aagaaa",
	"aatata",
	"aataca",
	"cataaa",
	"gataaa",
	"aatgaa",
	"tttaaa",
	"actaaa",
	"aataga",
	"aaaaag",
	"aaaaca",
	"ggggct"
};


/**
 * Vector with reversed hexamers. Same as hexamers, but reversed (not complemented). Used for searching with reverse_iterator.
 */
const std::vector<std::string> Utr3Finder::rHexamers = {
	"aaataa",
	"aaatta",
	"aaatat",
	"aaatga",
	"aaagaa",
	"atataa",
	"acataa",
	"aaatac",
	"aaatag",
	"aagtaa",
	"aaattt",
	"aaatca",
	"agataa",
	"gaaaaa",
	"acaaaa",
	"tcgggg"
};


/**
 * Constructor.
 */
Utr3Finder::Utr3Finder(const SeqStruct & tM) :
	txMut(tM)
{
	this->findPolyaMotif();
}


/**
 * Destructor.
 */
Utr3Finder::~Utr3Finder() {}


/**
 * Predicts the location of the UTR3's Poly(A) motif.
 */
void Utr3Finder::findPolyaMotif() {
	//transcript length is not always the same as sequence size
	size_t seqLength = txMut.txLength;
	size_t utr3Length = seqLength - this->txMut.utr3Start;
	size_t hitPos = Utr3Finder::noHitPos;
	//length of the sequence searched from transcript end
	size_t searchRange = 50;
	//ignoring nucleotides that aren't part of the transcript (possibly Poly(A) tail)
	size_t searchOffset = 0;
	if (txMut.seq.size() > seqLength) searchOffset = txMut.seq.size() - seqLength;

	BOOST_FOREACH (const std::string & pattern, this->rHexamers) {
		// The range in which the pattern is searched from the end of the transcript
		if (searchRange > utr3Length) searchRange = utr3Length; 

		auto posIt = std::search(this->txMut.seq.rbegin() + searchOffset, 
			this->txMut.seq.rbegin() + searchOffset + searchRange, 
			pattern.begin(), pattern.end());
		
		if (posIt != this->txMut.seq.rbegin() + searchOffset + searchRange) {
			/* Since positions are indexed starting with zero we have to substract 1 
			from the seqLength and the std::distance (+5 instead of +6) value */
			hitPos = (seqLength - 1) - (std::distance(this->txMut.seq.rbegin() + searchOffset, posIt) + 5);
			this->polyaMotifPos = hitPos;
			break;
		} 	
	}
	if (hitPos == Utr3Finder::noHitPos) this->polyaMotifPos = hitPos;
}


/**
 * Checks if variant is in the Poly(A) motif.
 */
bool Utr3Finder::isMutationInMotif() const {
	if (this->polyaMotifPos == Utr3Finder::noHitPos) return false;
	
	int diff = this->txMut.mutation->getMutPosition() + static_cast<int>(this->txMut.utr3Start) -
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
size_t Utr3Finder::getPolyaMotifPos() const {
	return this->polyaMotifPos;
}


/**
 * Returns the motif hexamer.
 */
std::string Utr3Finder::getMotifSequence() const {
	if (this->polyaMotifPos == Utr3Finder::noHitPos) return std::string();

	auto motifStart = this->txMut.seq.begin() + this->polyaMotifPos;
	auto motifEnd = this->txMut.seq.begin() + this->polyaMotifPos + 6;
	return std::string(motifStart, motifEnd);	
}


/**
 * Return string of the searched sequence.
 */
std::string Utr3Finder::getSequence() const {
	return this->txMut.seq;
}


/**
 * Return a string with information about the transcript and mutation.
 */
std::string Utr3Finder::writeInfo() const {
	std::stringstream ss;
	ss << *this->txMut.chrom 
		<< ", " << *this->txMut.genomicPos 
		<< ", " << *this->txMut.seqId
		<< ", " << *this->txMut.strand
		<< ", Poly(A) Pos: " << this->polyaMotifPos 
		<< ", utr3Start: " << this->txMut.utr3Start
		<< ", MutPos: " << this->txMut.utr3Start + txMut.mutation->getMutPosition()
		<< ", txLength: " << this->txMut.txLength;
	return ss.str();
}



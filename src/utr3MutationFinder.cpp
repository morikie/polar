#include <algorithm>
#include <sstream>
#include <vector>
#include <boost/foreach.hpp>
#include <seqan/find.h>
#include "utr3MutationFinder.hpp"


/**
 * Vector of known Poly(A) cleavage motifs to be not pathogenic. Sorted by frequency (highest first).
 */
const std::vector<std::string> Utr3MutationFinder::hexamers = {
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
	"aataac"
};


/**
 * Vector with reversed hexamers. Same as hexamers, just reversed.
 */
const std::vector<std::string> Utr3MutationFinder::rHexamers = {
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
	"caataa"
};


/**
 * Constructor.
 */
Utr3MutationFinder::Utr3MutationFinder(const TranscriptMutation & tM) :
	txMut(tM)
{
	this->findPolyaMotif();
}


/**
 * Destructor.
 */
Utr3MutationFinder::~Utr3MutationFinder() {}


/**
 * Predicts the location of the UTR3's Poly(A) motif.
 */
void Utr3MutationFinder::findPolyaMotif() {
	size_t seqLength = txMut.txLength;
	size_t hitPos = Utr3MutationFinder::noHitPos;
	//length of the sequence searched from transcript end
	size_t offset = 40;
	//ignoring the Poly(a) tail
	size_t searchOffset = 0;
	if (txMut.seq.size() > seqLength) searchOffset = txMut.seq.size() - seqLength;

	BOOST_FOREACH (const std::string & pattern, this->rHexamers) {
		size_t utr3Length = seqLength - this->txMut.utr3Start;
		// The range in which the pattern is searched from the end of the transcript
		if (offset > utr3Length) offset = utr3Length; 

		auto posIt = std::search(this->txMut.seq.rbegin() + searchOffset, this->txMut.seq.rbegin() + (searchOffset + offset), 
			pattern.begin(), pattern.end());
		
		if (posIt != this->txMut.seq.rbegin() + (searchOffset + offset)) {
			/* Since positions are indexed starting with zero we have to substract 1 
			from the seqLength and the std::distance (+5 instead of +6) value */
			hitPos = (seqLength - 1) - (std::distance(this->txMut.seq.rbegin() + searchOffset, posIt) + 5);
			this->polyaMotifPos = hitPos;
			break;
		} 	
	}
	if (hitPos == Utr3MutationFinder::noHitPos) this->polyaMotifPos = hitPos;
}


/**
 * Checks if variant is in the Poly(A) motif.
 */
bool Utr3MutationFinder::isMutationInMotif() const {
	if (this->polyaMotifPos == Utr3MutationFinder::noHitPos) return false;
	int diff = this->txMut.mutation.getMutPosition() + static_cast<int>(this->txMut.utr3Start) -
		static_cast<int>(this->polyaMotifPos);
//	std::cerr << "polyaMotifPos: " << this->polyaMotifPos << std::endl;
//	std::cerr << "strand: " << this->txMut.strand << std::endl;
//	std::cerr << "txName: " << this->txMut.seqId << std::endl;
//	std::cerr << "diff: " << diff << std::endl;
//	std::cerr << "txLength: " << this->txMut.seq.size() << std::endl;
//	std::cerr << "utr3Start: " << this->txMut.utr3Start << std::endl;
//	std::cerr << "mutationPos: " << this->txMut.mutation.getMutPosition() << std::endl; 
	if (diff < 6 && diff >= 0) {
		return true;
	} else {
		return false;
	}
}


/**
 * Returns the Poly(A) motif position in the transcript.
 */
size_t Utr3MutationFinder::getPolyaMotifPos() const {
	return this->polyaMotifPos;
}


/**
 * Returns the motif hexamer.
 */
std::string Utr3MutationFinder::getMotifSequence() const {
	if (this->polyaMotifPos == Utr3MutationFinder::noHitPos) return std::string();

	auto motifStart = this->txMut.seq.begin() + this->polyaMotifPos;
	auto motifEnd = this->txMut.seq.begin() + this->polyaMotifPos + 6;
	return std::string(motifStart, motifEnd);	
}


/**
 * Return a string with information about the transcript and mutation.
 */
std::string Utr3MutationFinder::writeLocation() const {
	std::stringstream ss;
	ss << *this->txMut.chrom 
		<< ", " << this->txMut.genomicPos 
		<< ", " << this->txMut.seqId 
		<< ", Poly(A) Pos: " << this->polyaMotifPos 
		<< ", utr3Start: " << this->txMut.utr3Start
		<< ", MutPos: " << this->txMut.utr3Start + txMut.mutation.getMutPosition();
	return ss.str();
}


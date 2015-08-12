#include <algorithm>
#include <sstream>
#include <vector>
#include <boost/foreach.hpp>
#include <seqan/find.h>
#include "utr3MutationFinder.hpp"


/**
 * Vector of known Poly(A) cleavage motifs to be not pathogenic.
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
	"aataga"
};


/**
 * Vector with reversed hexamers.
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
	"agataa"
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
 * Predicts the location of the UTR3's Poly(A) cleavage recognition motif.
 */
void Utr3MutationFinder::findPolyaMotif() {
	size_t seqLength = this->txMut.seq.size();
	size_t hitPos = Utr3MutationFinder::noHitPos;
	size_t offset = 40;	

	BOOST_FOREACH (const std::string & pattern, this->rHexamers) {
		size_t utr3Length = seqLength - this->txMut.utr3Start;
		// The range in which the pattern is searched from the end of the transcript
		if (offset > utr3Length) offset = utr3Length; 

		auto posIt = std::search(this->txMut.seq.rbegin(), this->txMut.seq.rbegin() + offset, 
			pattern.begin(), pattern.end());
		
		if (posIt != this->txMut.seq.rbegin() + offset) {
			/* Since positions are usually starting at zero we have to substract 1 
			from the seqLength and the std::distance (+5 instead of +6) value */
			hitPos = (seqLength - 1) - (std::distance(this->txMut.seq.rbegin(), posIt) + 5);
			this->polyaMotifPos = hitPos;
			break;
		} 	
	}
	if (hitPos == Utr3MutationFinder::noHitPos) this->polyaMotifPos = hitPos;
}


/**
 *
 */
bool Utr3MutationFinder::isMutationInMotif() const {
	if (this->polyaMotifPos == Utr3MutationFinder::noHitPos) return false;
	
	int diff = -1;
	diff = 	(this->txMut.mutation.getMutPosition() - 1) + static_cast<int>(this->txMut.utr3Start) -
		static_cast<int>(this->polyaMotifPos);
//	std::cerr << "polyaMotifPos: " << this->polyaMotifPos << std::endl;
//	std::cerr << "strand: " << this->txMut.strand << std::endl;
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
 *
 */
size_t Utr3MutationFinder::getPolyaMotifPos() const {
	return this->polyaMotifPos;
}


/**
 *
 */
std::string Utr3MutationFinder::getMotifSequence() const {
	if (this->polyaMotifPos == Utr3MutationFinder::noHitPos) return std::string();

	auto motifStart = this->txMut.seq.begin() + this->polyaMotifPos;
	auto motifEnd = this->txMut.seq.begin() + this->polyaMotifPos + 6;
	return std::string(motifStart, motifEnd);	
}

/**
 *
 */
std::string Utr3MutationFinder::writeLocation() const {
	std::stringstream ss;
	ss << this->txMut.chrom << ", " << this->txMut.genomicPos << ", " << this->txMut.seqId << ", Poly(A) Pos: " << this->polyaMotifPos << ", utr3Start: " << this->txMut.utr3Start;
	return ss.str();
}



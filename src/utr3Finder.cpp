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
	//cutoff is difference between true length (like str.size())
	//of the sequence and the designated length from external sources.
	size_t cutoff = txMut.seq.size() - seqLength;
	size_t utr3Length = seqLength - this->txMut.utr3Start;
	size_t hitPos = Utr3Finder::noHitPos;
	//length of the sequence searched from transcript end
	size_t searchRange = 100;
	//excepting remnants of the Poly(A) tail from the search
	size_t searchOffset = 0;
	auto it = this->txMut.seq.rbegin() + cutoff; 
	while (*it == 'a') {
		searchOffset++;
		it++;
	}
	if (searchRange > utr3Length) searchRange = utr3Length; 
	
	BOOST_FOREACH (const std::string & pattern, this->rHexamers) {
		//std::search returns an iterator to the start of the match
		//(in this case to the "end" of the match since we use reverse_iterators)
		auto posIt = std::search(this->txMut.seq.rbegin() + cutoff + searchOffset, 
			this->txMut.seq.rbegin() + cutoff + searchOffset + searchRange, 
			pattern.begin(), pattern.end());
		
		if (posIt != this->txMut.seq.rbegin() + cutoff + searchOffset + searchRange) {
			//Since positions are indexed starting with zero we have to substract 1 
			//from the seqLength and the std::distance (+5 instead of +6) value
			hitPos = (seqLength - 1) - (std::distance(this->txMut.seq.rbegin() + cutoff, posIt + 5));
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
	if (! this->txMut.mutation) return false;

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
 * Writes a string to the error console with information about the transcript and mutation (if available).
 */
void Utr3Finder::writeInfo() const {
	std::stringstream ss;
	ss << "Poly(A) Pos: " << this->polyaMotifPos 
		<< ", utr3Start: " << this->txMut.utr3Start
		<< ", txLength: " << this->txMut.txLength;
	if (this->txMut.chrom) ss << ", chrom: " <<  *this->txMut.chrom;
	if (this->txMut.genomicPos) ss << ", gPos: " <<  *this->txMut.genomicPos;
	if (this->txMut.seqId) ss << ", seqID: " << *this->txMut.seqId;
	if (this->txMut.strand) ss << ", strand: " << *this->txMut.strand;
	if (this->txMut.mutation) ss << ", MutPos: " << this->txMut.utr3Start + this->txMut.mutation->getMutPosition();
	
	std::cerr << ss.str() << std::endl;
}



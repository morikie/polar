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
	if (this->txMut.strand == "+") {	
		BOOST_FOREACH (const std::string & pattern, this->rHexamers) {
			size_t hitPos = 0;
			// The range in which the pattern is searched from the end of the transcript
			size_t offset = 100; //seqLength - this->txMut.utr3Start;
			if (offset > seqLength - this->txMut.utr3Start) offset = seqLength - this->txMut.utr3Start;

			auto posIt = std::search(this->txMut.seq.rbegin(), this->txMut.seq.rbegin() + offset, 
				pattern.begin(), pattern.end());

			if (posIt != this->txMut.seq.rbegin() + offset) {
				hitPos = seqLength - (std::distance(this->txMut.seq.rbegin(), posIt) + 6) ;
				this->polyaMotifPos = hitPos;
				break;
			}
			
		}
	} else {
		BOOST_FOREACH (const std::string & pattern, this->rHexamers) {
			size_t hitPos = 0;
			// The range in which the pattern is searched from the end of the transcript
			
			size_t offset = 100;
			if (offset > this->txMut.utr3Start) offset = this->txMut.utr3Start;

			auto posIt = std::search(this->txMut.seq.begin(), this->txMut.seq.begin() + offset, 
				pattern.begin(), pattern.end());
			
			if (posIt != this->txMut.seq.begin() + offset) {
				hitPos = std::distance(this->txMut.seq.begin(), posIt ) + 6;
				if (hitPos != this->txMut.utr3Start + 6) {
					this->polyaMotifPos = hitPos;
					break;
				}
			}
		}
	}
}


/**
 *
 */
bool Utr3MutationFinder::isMutationInMotif() const {
	if (this->polyaMotifPos == 0) return false;
	
	int diff = -1;
	if (this->txMut.strand == "+") {
		std::cerr << "in '+' strand" << std::endl;
		diff = static_cast<int>(this->polyaMotifPos) - 
		(static_cast<int>(this->txMut.mutation.utr3MutPos) + static_cast<int>(this->txMut.utr3Start));
	} else {
		diff = (static_cast<int>(this->polyaMotifPos) - 
		(static_cast<int>(this->txMut.utr3Start)) + static_cast<int>(this->txMut.mutation.utr3MutPos));
	}
	
	std::cerr << "polyaMotifPos: " << this->polyaMotifPos << std::endl;
	std::cerr << "strand: " << this->txMut.strand << std::endl;
	std::cerr << "diff: " << diff << std::endl;
	std::cerr << "txLength: " << this->txMut.seq.size() << std::endl;
	std::cerr << "utr3Start: " << this->txMut.utr3Start << std::endl;
	std::cerr << "mutationPos: " << this->txMut.mutation.utr3MutPos << std::endl; 
	if (diff < 6 && diff >= 0) {
		return true;
	} else {
		return false;
	}
}


/**
 *
 */
std::string Utr3MutationFinder::writeLocation() const {
	std::stringstream ss;
	ss << this->txMut.chrom << ", " << this->txMut.genomicPos << ", " << this->txMut.seqId;
	return ss.str();
}


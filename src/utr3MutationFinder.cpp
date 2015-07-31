#include <sstream>
#include <vector>
#include <boost/foreach.hpp>
#include <seqan/find.h>
#include "utr3MutationFinder.hpp"


/**
 * Set of known Poly(A) cleavage motifs to be not pathogenic.
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
 * Constructor,
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
//	seqan::Finder<const std::string> finder(this->sequence);
//	
//	std::string p = "AATAAA";
//	seqan::Pattern<std::string, seqan::Pex<seqan::Hierarchical, seqan::AhoCorasick> > pat(p, -2);
//	std::cerr << "before find(...)" << std::endl;
//	while (seqan::find(finder, pat)) {
//		while (findBegin(finder, pat)){
//			std::cerr << infix(finder) << std::endl;
//		}
//	}
	
	size_t seqLength = this->txMut.seq.size();
	if (seqLength > 60) {	
//		std::cerr << this->txMut.seqId << ":" << std::endl;
//		std::cerr << "this->txMut.utr3Start: " << this->txMut.utr3Start << std::endl;
//		std::cerr << "seqLength: " << seqLength << std::endl;
//		std::cerr << "this->txMut.seq: " << this->txMut.seq << std::endl;
		
		BOOST_FOREACH (const std::string & pattern, this->hexamers) {
			size_t hitPos = txMut.utr3Start;	
			while ((hitPos = this->txMut.seq.find(pattern, hitPos + 1)) != std::string::npos) {
				this->polyaMotifPos = hitPos;
				std::cerr << "Found hit of " << pattern << " at " << hitPos << std::endl;	
			}
			break;
		}
	} else {
		std::cerr << this->txMut.seq << std::endl;
	}
}


/**
 *
 */
bool Utr3MutationFinder::isMutationInMotif() const {
	int diff = static_cast<int>(this->polyaMotifPos) - (static_cast<int>(this->txMut.mutation.utr3MutPos) + static_cast<int>(this->txMut.utr3Start));
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


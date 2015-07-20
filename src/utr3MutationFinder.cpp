#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#include <seqan/find.h>
#include "utr3MutationFinder.hpp"


const std::set<std::string> MutationFinder::hexamers = boost::assign::list_of
	("AATAAA")
	("ATTAAA")
	("TATAAA")
	("AGTAAA")
	("AAGAAA")
	("AATATA")
	("AATACA")
	("CATAAA")
	("GATAAA")
	("AATGAA")
	("TTTAAA")
	("ACTAAA")
	("AATAGA")
;


MutationFinder::MutationFinder(const std::string & seq) :
sequence(seq)
{
	this->findConsensus();
}


MutationFinder::~MutationFinder() {}


void MutationFinder::findConsensus() {
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
	
	size_t seqLength = this->sequence.size();
	if (seqLength > 40) {
		BOOST_FOREACH (std::string pattern, this->hexamers) {
			size_t hitPos = 0;
			if ((hitPos = this->sequence.find(pattern, seqLength - 40)) != std::string::npos) {
				std::cerr << "Found hit of " << pattern << " at " << hitPos << std::endl;
			}
			
		}
	} else {
		std::cerr << this->sequence << std::endl;
	}

	
	

}

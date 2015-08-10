#include <iostream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/filesystem/path.hpp>
#include "knownGeneParser.hpp"
#include "jannovarVcfParser.hpp"
#include "readTranscripts.hpp"
#include "transcriptMutation.hpp"
#include "readTranscriptMutation.hpp"
#include "utr3MutationFinder.hpp"
#include "polar.hpp"

namespace fs = boost::filesystem;

int main (int args, char * argv[]) {
	fs::path knownGene = "knownGene.txt";
	fs::path transcripts = "knownGeneMrna.txt";
	fs::path vcfFile = "vcf-example.jv.vcf";
	KnownGeneParser txValues(knownGene);
	ReadTranscripts tx(transcripts);
	JannovarVcfParser variants(vcfFile);
	std::vector<Utr3MutationFinder> utr3MutFinderVector;	
	std::vector<TranscriptMutation> transMutVector;
	
	readTranscriptMutation (transMutVector, variants, txValues, tx);
	
	try {
		for (auto it = transMutVector.begin(); it != transMutVector.end(); ++it) {
			utr3MutFinderVector.push_back(Utr3MutationFinder(*it));
		}
	} catch (...) {
		std::cerr << "error occured" << std::endl;
	}
	
	//Amount of UTR3 mutations which do not affect the Poly(A) motif
	size_t noFinds = 0;
	BOOST_FOREACH (const Utr3MutationFinder & utr3MutFi, utr3MutFinderVector) {
		if (utr3MutFi.isMutationInMotif()) {
			std::cerr << "Poly(A) motif mutation detected: " << utr3MutFi.writeLocation() << std::endl;
		} else {
			noFinds++;
		}
	}

	size_t undetectedUtr3Motifs = 0;
	
	BOOST_FOREACH (const Utr3MutationFinder & utr3MutFi, utr3MutFinderVector) {
		if (utr3MutFi.getPolyaMotifPos() == Utr3MutationFinder::noHitPos) {
			undetectedUtr3Motifs++;
		}
	}
	std::cerr << "undetectedUtr3Motifs: " << undetectedUtr3Motifs << std::endl;
	std::cerr << "utr3MutFinderVector.size(): " << utr3MutFinderVector.size() << std::endl;
	std::cerr << "noFinds: " << noFinds << std::endl;

	return EXIT_SUCCESS;
}

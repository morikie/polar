#include <iostream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
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

	std::vector<TranscriptMutation> transMutVector;

	readTranscriptMutation(transMutVector, variants, txValues, tx);

	std::vector<Utr3MutationFinder> utr3MutFinderVector;	
	try {

		for (auto it = transMutVector.begin(); it != transMutVector.end(); ++it) {
			utr3MutFinderVector.push_back(Utr3MutationFinder(*it));
		}
	} catch (...) {
		std::cerr << "error occured" << std::endl;
	}

	BOOST_FOREACH(const Utr3MutationFinder & utr3MutFi, utr3MutFinderVector) {
		if (utr3MutFi.isMutationInMotif()) {
			std::cerr << "Poly(A) motif mutation detected!" << std::endl;
		} else {
			std::cerr << "Mutation doesn't affect Poly(A) site." << std::endl;
		}
	}

	return EXIT_SUCCESS;
}

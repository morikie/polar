#include <iostream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/filesystem/path.hpp>
#include "knownGeneParser.hpp"
#include "jannovarVcfParser.hpp"
#include "readTranscripts.hpp"
#include "seqStruct.hpp"
#include "readSeqStruct.hpp"
#include "utr3Finder.hpp"
#include "polar.hpp"

namespace fs = boost::filesystem;

int main (int args, char * argv[]) {
	fs::path knownGene = "knownGene.txt";
	fs::path transcripts = "knownGeneMrna.txt";
	fs::path vcfFile = "vcf-example.jv.vcf";
	KnownGeneParser txValues(knownGene);
	ReadTranscripts tx(transcripts);
	JannovarVcfParser variants(vcfFile);
	std::vector<Utr3Finder> utr3MutFinderVector;	
	std::vector<SeqStruct> transMutVector;
	
	readSeqStruct (transMutVector, variants, txValues, tx);
	
	try {
		for (auto it = transMutVector.begin(); it != transMutVector.end(); ++it) {
			utr3MutFinderVector.push_back(Utr3Finder(*it));
		}
	} catch (...) {
		std::cerr << "error occured" << std::endl;
	}
	
	//Amount of UTR3 mutations which do not affect the Poly(A) motif
	size_t noFinds = 0;
	BOOST_FOREACH (const Utr3Finder & utr3MutFi, utr3MutFinderVector) {
		if (utr3MutFi.isMutationInMotif()) {
			std::cerr << "Poly(A) motif mutation detected: " << utr3MutFi.writeInfo() << std::endl;
		} else {
			noFinds++;
		}
	}

	size_t undetectedUtr3Motifs = 0;

	BOOST_FOREACH (const Utr3Finder & utr3MutFi, utr3MutFinderVector) {
		if (utr3MutFi.getPolyaMotifPos() == Utr3Finder::noHitPos) {
			std::cerr << "couldn't find motif: " << utr3MutFi.writeInfo() << std::endl;
			undetectedUtr3Motifs++;
		}
	}

	std::cerr << "undetectedUtr3Motifs: " << undetectedUtr3Motifs << std::endl;
	std::cerr << "utr3MutFinderVector.size(): " << utr3MutFinderVector.size() << std::endl;
	std::cerr << "noFinds: " << noFinds << std::endl;

	return EXIT_SUCCESS;
}

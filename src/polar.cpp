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
#include "utr3FinderNaive.hpp"
#include "utr3FinderFuzzy.hpp"
#include "polar.hpp"

namespace fs = boost::filesystem;

int main (int args, char * argv[]) {
	fs::path knownGene = "ucsc_data/knownGene.txt";
	fs::path transcripts = "ucsc_data/knownGeneTxMrna.txt";
	//fs::path vcfFile = "vcf/vcf-example.jv.vcf";
	fs::path vcfFile = "vcf/1000Genomes/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.jv.vcf";
	//fs::path vcfFile = "vcf/1000Genomes/taddaa2.txt";
	//fs::path vcfFile = "vcf/1000Genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.jv.vcf";
	//fs::path vcfFile = "vcf/vcf-example.jv.vcf";
	KnownGeneParser txValues(knownGene);
	ReadTranscripts tx(transcripts);
	JannovarVcfParser variants(vcfFile);
	std::vector<Utr3Finder*> utr3MutFinderVector;	
	std::vector<SeqStruct> transMutVector;
	
	readSeqStruct (transMutVector, variants, txValues, tx);
	

	try {
		for (auto it = transMutVector.begin(); it != transMutVector.end(); ++it) {
			utr3MutFinderVector.push_back(new Utr3FinderNaive(*it));
		}
	} catch (...) {
	//	std::cerr << "error occured" << std::endl;
	}
	
	//Amount of UTR3 mutations which do not affect the Poly(A) motif
	size_t noPolya = 0;
	BOOST_FOREACH (const Utr3Finder * utr3MutFi, utr3MutFinderVector) {
		if (utr3MutFi->isMutationInMotif()) {
			std::cerr << "Poly(A) motif mutation detected: ";
			utr3MutFi->writeInfo();
		} else {
			noPolya++;
		}
	}

	size_t undetectedUtr3Motifs = 0;
	BOOST_FOREACH (const Utr3Finder * utr3MutFi, utr3MutFinderVector) {
		if (utr3MutFi->getPolyaMotifPos()[0] == Utr3Finder::noHitPos) {
			//tr3MutFi.writeInfo();
			//std::cerr << utr3MutFi.getSequence() << std::endl;
			//std::cerr << "couldn't find motif: " << utr3MutFi.writeInfo() << std::endl;
			undetectedUtr3Motifs++;
		}
	}

	std::cerr << "undetectedUtr3Motifs: " << undetectedUtr3Motifs << std::endl;
	std::cerr << "utr3MutFinderVector.size(): " << utr3MutFinderVector.size() << std::endl;
	std::cerr << "noPolya: " << noPolya << std::endl;

	return EXIT_SUCCESS;
}


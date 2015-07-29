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

	std::vector<TranscriptMutation> transMutVector;

	readTranscriptMutation(transMutVector, variants, txValues, tx);
	
	std::cerr << "transMutVector.size(): " << transMutVector.size() << std::endl;
	std::cerr << transMutVector[0].mutation.utr3MutPos << std::endl;

	return EXIT_SUCCESS;
}

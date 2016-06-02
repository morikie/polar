#include <iostream>
#include <string>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <seqan/seq_io.h>

#include "buildIndexFile.hpp"
#include "knownGeneParser.hpp"
#include "jannovarVcfParser.hpp"
#include "readTranscripts.hpp"
#include "seqStruct.hpp"
#include "buildSeqStruct.hpp"
#include "utr3Finder.hpp"
#include "utr3FinderNaive.hpp"
#include "utr3FinderFuzzy.hpp"
#include "polar.hpp"

extern "C" {
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/utils.h>
}


namespace fs = boost::filesystem;

int main (int args, char * argv[]) {
	const char * test = "cgtgtgattgcacg";
	char *mfe_structure = static_cast<char*>(vrna_alloc(sizeof(char) * (strlen(test) + 1)));
	vrna_fold_compound_t * test_fold = vrna_fold_compound (test, NULL, VRNA_OPTION_PF | VRNA_OPTION_MFE);
	vrna_mfe(test_fold, mfe_structure);
	std::cerr << test << std::endl << mfe_structure << std::endl;
	vrna_pf(test_fold, NULL);
	std::cerr << "exp_matrices->length: " <<  test_fold->exp_matrices->length << std::endl;
	vrna_plist_t * pl = vrna_plist_from_probs(test_fold, 1e-5);

	for (size_t k = 0; pl[k].i>0 || pl[k].j>0; k++) {
		std::cerr << pl[k].i << ", " << pl[k].j << ": " << pl[k].p << std::endl;
	}
	std::cerr << __FUNCTION__ << std::endl;
	fs::path knownGene = "ucsc_data/knownGene.txt";
	fs::path transcripts = "ucsc_data/knownGeneTxMrna.txt";
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	//fs::path vcfFile = "vcf/sensitivity.jv.vcf";
	fs::path vcfFile = "vcf/vcf-example.jv.vcf";
	//fs::path vcfFile = "vcf/1000Genomes/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.jv.vcf";
	//fs::path vcfFile = "vcf/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.jv.vcf";
	//fs::path vcfFile = "vcf/1000Genomes/taddaa2.txt";
	//fs::path vcfFile = "vcf/1000Genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.jv.vcf";
	KnownGeneParser txValues(knownGene);
	ReadTranscripts tx(transcripts);
	JannovarVcfParser variants(vcfFile);
	std::vector<std::shared_ptr<Utr3Finder> > utr3FinderVector;	
	std::vector<SeqStruct> seqStructVector;
	
	
	if (args > 1) {
		std::cerr << "Naive" << std::endl;
		buildSeqStructFromTranscripts (seqStructVector, variants, txValues, tx);
		std::cerr << seqStructVector.size() << std::endl;
		try {
			for (auto it = seqStructVector.begin(); it != seqStructVector.end(); ++it) {
				utr3FinderVector.push_back(std::shared_ptr<Utr3Finder>(new Utr3FinderNaive(*it)));
			}
		} catch (...) {
			std::cerr << "error occured creating Utr3Finder object" << std::endl;
		}
		
		//Amount of UTR3 mutations which do not affect the Poly(A) motif
		size_t noPolya = 0;
		BOOST_FOREACH (std::shared_ptr<Utr3Finder> utr3MutFi, utr3FinderVector) {
			if (utr3MutFi->isMutationInMotif()) {
				std::cerr << "Poly(A) motif mutation detected: ";
				utr3MutFi->writeInfo();
			} else {
				noPolya++;
			}
		}

		size_t undetectedUtr3Motifs = 0;
		BOOST_FOREACH (std::shared_ptr<Utr3Finder> utr3MutFi, utr3FinderVector) {
			if (utr3MutFi->getPolyaMotifPos()[0].pos == Utr3Finder::noHitPos) {
				//tr3MutFi.writeInfo();
				//std::cerr << utr3MutFi.getSequence() << std::endl;
				//std::cerr << "couldn't find motif: " << utr3MutFi.writeInfo() << std::endl;
				undetectedUtr3Motifs++;
			}
		}

		std::cerr << "undetectedUtr3Motifs: " << undetectedUtr3Motifs << std::endl;
		std::cerr << "utr3FinderVector.size(): " << utr3FinderVector.size() << std::endl;
		std::cerr << "noPolya: " << noPolya << std::endl;

	} else {
		std::cerr << "Fuzzy" << std::endl;
		if (! fs::exists(refGenomeIndex)) {
			buildIndexFile(referenceGenome);	
		}
		
		seqan::FaiIndex faiIndex;
		if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
			std::cerr << "could not open index file for " << referenceGenome << std::endl;
		}

		buildSeqStructFromGenome(seqStructVector, variants, faiIndex); 
		
		try {
			for (auto it = seqStructVector.begin(); it != seqStructVector.end(); ++it) {
				utr3FinderVector.push_back(std::shared_ptr<Utr3Finder>(new Utr3FinderFuzzy(*it)));
			}
		} catch (...) {
			std::cerr << "error occured creating Utr3Finder object" << std::endl;
		}

		BOOST_FOREACH (std::shared_ptr<Utr3Finder> utr3Fuzzy, utr3FinderVector) {
			size_t seqLength = utr3Fuzzy->getSequence().size();
			if (! utr3Fuzzy->getPolyaMotifPos().empty()) { 
				std::vector<Utr3Finder::Utr3FinderResult> positions = utr3Fuzzy->getPolyaMotifPos();
				BOOST_FOREACH(const Utr3Finder::Utr3FinderResult & result, positions) {
					if (result.strand == "+" && result.pos >= 94 && result.pos <= 99) {
						utr3Fuzzy->writeInfo();
					} else if (result.strand == "-" && result.pos >= 99 && result.pos <= 104) {
						utr3Fuzzy->writeInfo();	
					}
				}
			}	
		}
	}

	return EXIT_SUCCESS;
}


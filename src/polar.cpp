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

namespace fs = boost::filesystem;


int main (int args, char * argv[]) {
	std::cerr << __FUNCTION__ << std::endl;
	fs::path knownGene = "ucsc_data/knownGene.txt";
	fs::path transcripts = "ucsc_data/knownGeneTxMrna.txt";
	fs::path vcfFile = "vcf/vcf-example.jv.vcf";
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	//fs::path vcfFile = "vcf/1000Genomes/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.jv.vcf";
	//fs::path vcfFile = "vcf/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.jv.vcf";
	//fs::path vcfFile = "vcf/1000Genomes/taddaa2.txt";
	//fs::path vcfFile = "vcf/1000Genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.jv.vcf";
	//fs::path vcfFile = "vcf/vcf-example.jv.vcf";
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
			if (utr3MutFi->getPolyaMotifPos()[0] == Utr3Finder::noHitPos) {
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
		
		seqan::CharString test;
		seqan::readRegion(test, faiIndex, 22, 49106919 - 150, 49106919 + 150);
		std::cerr << test << std::endl;

		buildSeqStructFromGenome(seqStructVector, variants, faiIndex); 
		
		try {
			for (auto it = seqStructVector.begin(); it != seqStructVector.end(); ++it) {
				utr3FinderVector.push_back(std::shared_ptr<Utr3Finder>(new Utr3FinderFuzzy(*it)));
			}
		} catch (...) {
			std::cerr << "error occured creating Utr3Finder object" << std::endl;
		}

		BOOST_FOREACH (std::shared_ptr<Utr3Finder> utr3Fuzzy, utr3FinderVector) {
				if (! utr3Fuzzy->getPolyaMotifPos().empty()) { 
					std::vector<size_t> positions = utr3Fuzzy->getPolyaMotifPos();
					BOOST_FOREACH(size_t pos, positions) {
						if (pos >= 144 && pos <= 150) {
							utr3Fuzzy->writeInfo();
						}
					}
				}
		}
	}

	return EXIT_SUCCESS;
}


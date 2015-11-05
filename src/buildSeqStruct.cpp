#include <regex>
#include <sstream>
#include <boost/optional.hpp>
#include <seqan/seq_io.h>
#include "buildSeqStruct.hpp"
#include "hgvsParser.hpp"


/**
 * Calculate the 3' UTR start. Keep in mind, that the start is not always in the last exon (i.e. the UTR might span over several exons).
 */
inline size_t findUtrStart (const TxProperties & txProp, const size_t & txLength) {
	size_t utr3Length = 0;
	if (txProp.strand == "+") {

		for (int i = txProp.exonStarts.size() - 1; i >= 0; i--) {
			if (txProp.cdsEnd >= txProp.exonStarts[i]) {
				utr3Length = txProp.exonEnds[i] - txProp.cdsEnd;
				for (size_t j = static_cast<size_t>(i) + 1; j < txProp.exonStarts.size(); j++) { 
					utr3Length += txProp.exonEnds[j] - txProp.exonStarts[j];
				}
				break;
			}
		}
	} else {
		for (size_t i = 0; i < txProp.exonEnds.size(); i++) {
			if (txProp.cdsStart <= txProp.exonEnds[i]) {
				utr3Length = txProp.cdsStart - txProp.exonStarts[i];
				for (int j = static_cast<int>(i) - 1; j >= 0; j--) { 
					utr3Length += txProp.exonEnds[j] - txProp.exonStarts[j];
				}
				break;
			}
		}
	}
	return txLength - utr3Length;
}


/**
 * Calculate length of transcript by summing over the exon lengths.
 */
inline size_t getTxLength (const TxProperties & txProp) {
	if (txProp.exonStarts.size() != txProp.exonStarts.size()) return 0;
	
	size_t sum = 0;
	for (size_t i = 0; i < txProp.exonStarts.size(); i++) {
		sum += txProp.exonEnds[i] - txProp.exonStarts[i];
	}
	return sum;

}


/**
 * Fills up the SeqStruct objects needed for the Utr3Finder.
 */
bool buildSeqStructFromTranscripts (std::vector<SeqStruct> & transMutVector, 
	const JannovarVcfParser & vcfParser, 
	const KnownGeneParser & txValues, 
	const ReadTranscripts & txSequences) 
{
	auto itBegin = (vcfParser.getData()).begin();
	auto itEnd = (vcfParser.getData()).end();
	
	size_t inIntron = 0;
	size_t count = 0;
	size_t count_affectedTx = 0;
	size_t count_affectedUtr3 = 0;
	//size_t counter = 1;
	for ( ; itBegin != itEnd; itBegin++) {
		count++;
		BOOST_FOREACH(const vcfTranscripts & vcfTx, itBegin->second) {
			count_affectedTx++;
			if (vcfTx.jvVariantType.find("3_prime_utr_variant") != std::string::npos) {	
				//std::stringstream ss;
				//ss << "3_prime_utr_variant" << counter;
				/*if (vcfTx.jvVariantType.find(ss.str()) != std::string::npos) {
					std::cerr << ss.str() << " wasn't found!" << std::endl;
				}
				counter++;*/
				count_affectedUtr3++;
				auto & mapValues = txValues.getValueByKey(vcfTx.txName);
				size_t txLength = getTxLength(mapValues); 
				auto & chrom = (itBegin->first).first;
				auto & genePos = (itBegin->first).second;
				auto & strand = mapValues.strand;

				//ignore transcripts that have no coding sequence 
				if (mapValues.cdsEnd == mapValues.cdsStart) continue;

				size_t utr3Start = findUtrStart(mapValues, txLength);
				SeqStruct transMut = {
					txSequences.getValueByKey(vcfTx.txName),
					utr3Start,
					txLength,
					boost::optional<const HgvsParser>(HgvsParser(vcfTx.hgvsString)),
					boost::optional<const std::string &>(chrom),
					boost::optional<const size_t &>(genePos),
					boost::optional<const std::string &>(strand),
					boost::optional<const std::string &>(vcfTx.txName)
				};
				if (! transMut.mutation->isIntronic()) {
					transMutVector.push_back(transMut);
				} else {
					inIntron++;		
				}
			}
		}
	}
//	std::cerr << "jv Tx counts: " << count_affectedTx << std::endl;
//	std::cerr << "jv UTR3 counts: " << count_affectedUtr3 << std::endl;
//	std::cerr << "UTR variants in introns: " << inIntron << std::endl;
	return true;
}


/**
 * Fills up the SeqStruct objects needed for the Utr3Finder.
 */
bool buildSeqStructFromGenome (std::vector<SeqStruct> & transMutVector, 
	const JannovarVcfParser & vcfParser,
	const seqan::FaiIndex & fai)
{
	auto itBegin = (vcfParser.getData()).begin();
	auto itEnd = (vcfParser.getData()).end();
	
	for ( ; itBegin != itEnd; itBegin++) {
		BOOST_FOREACH(const vcfTranscripts & vcfTx, itBegin->second) {
			if (vcfTx.jvVariantType.find("3_prime_utr_variant") != std::string::npos) {	
				auto & chrom = (itBegin->first).first;
				auto & genePos = (itBegin->first).second;
				
				std::string chr = chrom;
				chr.erase(0,3);
				unsigned idx;
				if (chr == "X") {
					idx = 22;
				} else if (chr == "Y") {
					idx = 23;
				} else {
					idx = std::stoi(chr);
				}
				

				seqan::String<char> seq;
				
				size_t startRange;
				size_t endRange;
				if (genePos < 150) {
					startRange = 0;
				} else {
					startRange = genePos - 150;
				}
				if (genePos + 150 > 

				seqan::readRegion(seq, fai, idx, startRange, genePos + 150);
				std::string sequence(seqan::toCString(seq));


				SeqStruct transMut = {
					sequence,
					boost::none,
					boost::none,
					boost::optional<const HgvsParser>(HgvsParser(vcfTx.hgvsString)),
					boost::optional<const std::string &>(chr),
					boost::optional<const size_t &>(genePos),
					boost::none,
					boost::optional<const std::string &>(vcfTx.txName)
				};
				if (! transMut.mutation->isIntronic()) 
					transMutVector.push_back(transMut);
			}
		}
	}
	return true;
}

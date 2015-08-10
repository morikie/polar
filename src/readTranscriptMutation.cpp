#include <numeric>
#include "readTranscriptMutation.hpp"
#include "hgvsParser.hpp"


inline size_t findUtrStart (const TxProperties & txProp, const size_t & txLength) {
	if (txProp.strand == "+") {
		size_t utr3Length = 0;

		for (int i = txProp.exonStarts.size() - 1; i >= 0; i--) {
			if (txProp.cdsEnd >= txProp.exonStarts[i]) {
				utr3Length = txProp.exonEnds[i] - txProp.cdsEnd;
				for (size_t j = static_cast<size_t>(i) + 1; j < txProp.exonStarts.size(); j++) { 
					utr3Length += txProp.exonEnds[j] - txProp.exonStarts[j];
				}
				break;
			}
		}
		return txLength - utr3Length;
	} else {

		size_t utr3Length = 0;

		for (size_t i = 0; i < txProp.exonEnds.size(); i++) {
			if (txProp.cdsStart <= txProp.exonEnds[i]) {
				utr3Length = txProp.cdsStart - txProp.exonStarts[i];
				for (int j = static_cast<int>(i) - 1; j >= 0; j--) { 
					utr3Length += txProp.exonEnds[j] - txProp.exonStarts[j];
				}
				break;
			}
		}
		return utr3Length;
	}
}


/**
 *
 */
bool readTranscriptMutation (std::vector<TranscriptMutation> & transMutVector, 
	const JannovarVcfParser & vcfParser, 
	const KnownGeneParser & txValues, 
	const ReadTranscripts & txSequences) 
{
	auto itBegin = (vcfParser.getData()).begin();
	auto itEnd = (vcfParser.getData()).end();
	
	size_t inIntron = 0;
	for ( ; itBegin != itEnd; itBegin++) {
			
		BOOST_FOREACH(const vcfTranscripts & vcfTx, itBegin->second) {
			if (vcfTx.jvVariantType == "3_prime_utr_variant") {
				auto & mapValues = txValues.getValueByKey(vcfTx.txName);
				size_t txLength = txSequences.getValueByKey(vcfTx.txName).size();
				auto & strand = txValues.getValueByKey(vcfTx.txName).strand;
				
				if (mapValues.cdsEnd == mapValues.cdsStart) continue;

				size_t utr3Start = findUtrStart(mapValues, txLength);

				TranscriptMutation transMut = {
					(itBegin->first).first,
					(itBegin->first).second,
					strand,
					vcfTx.txName,
					txSequences.getValueByKey(vcfTx.txName),
					HgvsParser(vcfTx.hgvsString),
					utr3Start
				};
				if (! transMut.mutation.isIntronic()) {
					transMutVector.push_back(transMut);
				} else {
					inIntron++;		
				}

			}
		}
	}
	std::cerr << "UTR variants in introns: " << inIntron << std::endl;
	return true;
}

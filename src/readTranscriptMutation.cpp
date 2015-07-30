#include "readTranscriptMutation.hpp"
#include "hgvsParser.hpp"

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

	for ( ; itBegin != itEnd; itBegin++) {
		
		BOOST_FOREACH(const vcfTranscripts & vcfTx, itBegin->second) {
			if (vcfTx.jvVariantType == "3_prime_utr_variant") {
				auto & mapValues = txValues.getValueByKey(vcfTx.txName);
				size_t txLength = txSequences.getValueByKey(vcfTx.txName).size();
				unsigned int utr3Start = txLength - (mapValues.txEnd - mapValues.cdsEnd);
				

				TranscriptMutation transMut = {
					vcfTx.txName,
					txSequences.getValueByKey(vcfTx.txName),
					HgvsParser(vcfTx.hgvsString),
					utr3Start
				};
				transMutVector.push_back(transMut);
			}
		}
	}
	return true;
}

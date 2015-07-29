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
		
		BOOST_FOREACH(vcfTranscripts vcfTx, itBegin->second) {
			if (vcfTx.jvVariantType == "3_prime_utr_variant") {
				
				TranscriptMutation transMut = {
					vcfTx.txName,
					txSequences.getValueByKey(vcfTx.txName),
					HgvsParser(vcfTx.hgvsString),
					txValues.getValueByKey(vcfTx.txName).cdsEnd
				};
				transMutVector.push_back(transMut);
			}
		}
	}
	return true;
}

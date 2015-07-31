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
				auto & strand = txValues.getValueByKey(vcfTx.txName).strand;
				
				unsigned int utr3Start = 0;
				if (strand == "+") { 
					utr3Start = txLength - (mapValues.txEnd - mapValues.cdsEnd);
				} else {
					utr3Start = mapValues.cdsStart - mapValues.txStart;
				}

				TranscriptMutation transMut = {
					(itBegin->first).first,
					(itBegin->first).second,
					strand,
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

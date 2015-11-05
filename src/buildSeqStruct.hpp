#ifndef __READTRANSCRIPTMUATION_HPP__
#define __READTRANSCRIPTMUATION_HPP__

#include <vector>
#include <seqan/seq_io.h>
#include "jannovarVcfParser.hpp"
#include "knownGeneParser.hpp"
#include "readTranscripts.hpp"
#include "seqStruct.hpp"


size_t findUtr3Start(const TxProperties & txProp, const size_t & txLength);


size_t getTxLength(const TxProperties & txProp);


/**
 *
 */
bool buildSeqStructFromTranscripts (std::vector<SeqStruct> & transMutVector, 
	const JannovarVcfParser & vcfParser, 
	const KnownGeneParser & txValues, 
	const ReadTranscripts & txSequences);


bool buildSeqStructFromGenome(std::vector<SeqStruct> & seqStructVector,
	const JannovarVcfParser & vcfParser,
	const seqan::FaiIndex & fai);


#endif /* __READTRANSCRIPTMUATION_HPP__ */


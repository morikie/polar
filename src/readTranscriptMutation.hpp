#ifndef __READTRANSCRIPTMUATION_HPP__
#define __READTRANSCRIPTMUATION_HPP__

#include <vector>
#include "jannovarVcfParser.hpp"
#include "knownGeneParser.hpp"
#include "readTranscripts.hpp"
#include "transcriptMutation.hpp"

/**
 *
 */
bool readTranscriptMutation (std::vector<TranscriptMutation> & transMutVector, 
	const JannovarVcfParser & vcfParser, 
	const KnownGeneParser & txValues, 
	const ReadTranscripts & txSequences);

#endif /* __READTRANSCRIPTMUATION_HPP__ */

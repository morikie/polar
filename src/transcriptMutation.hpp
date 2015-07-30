#ifndef __TRANSCRIPTMUTATION_HPP__
#define __TRANSCRIPTMUTATION_HPP__

#include <string>
#include "hgvsParser.hpp"

struct TranscriptMutation {
	
	const std::string & seqId;
	const std::string & seq;
	const HgvsParser mutation;

	size_t utr3Start;

};

#endif /* __TRANSCRIPTMUTATION_HPP__ */


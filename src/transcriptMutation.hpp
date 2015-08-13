#ifndef __TRANSCRIPTMUTATION_HPP__
#define __TRANSCRIPTMUTATION_HPP__

#include <string>
#include "hgvsParser.hpp"

struct TranscriptMutation {
	
	const std::string & chrom;
	const unsigned int & genomicPos;
	const std::string strand;

	const std::string & seqId;
	const std::string & seq;
	const HgvsParser mutation;

	const size_t utr3Start;
	const size_t txLength;
};

#endif /* __TRANSCRIPTMUTATION_HPP__ */


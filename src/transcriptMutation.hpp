#ifndef __TRANSCRIPTMUTATION_HPP__
#define __TRANSCRIPTMUTATION_HPP__

#include <string>
#include <boost/optional.hpp>
#include "hgvsParser.hpp"


struct TranscriptMutation {	
	const boost::optional<std::string> & chrom;
	const size_t & genomicPos;
	const std::string strand;

	const std::string & seqId;
	const std::string & seq;
	const HgvsParser mutation;

	const size_t utr3Start;
	const size_t txLength;
};

#endif /* __TRANSCRIPTMUTATION_HPP__ */


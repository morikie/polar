#ifndef __TRANSCRIPTMUTATION_HPP__
#define __TRANSCRIPTMUTATION_HPP__

#include <string>
#include <boost/optional.hpp>
#include "hgvsParser.hpp"


/**
 * This class is used as a wrapper for the arguments "needed" by Utr3Finder  
 */
struct SeqStruct {	
	const std::string & seq;
	boost::optional<const size_t> utr3Start;
	boost::optional<const size_t> txLength;
	//HGVS string parser
	boost::optional<const HgvsParser> mutation;
	
	//information about the sequence
	boost::optional<const std::string &> chrom;
	boost::optional<const size_t &> genomicPos;
	boost::optional<const std::string &> strand;
	
	//name of the sequence (e.g. transcript ID)
	boost::optional<const std::string &> seqId;
};

#endif /* __TRANSCRIPTMUTATION_HPP__ */


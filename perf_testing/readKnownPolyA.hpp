#ifndef __READKNOWNPOLYA_HPP__
#define __READKNOWNPOLYA_HPP__

#include <fstream>
#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;


struct KnownPolyA {
	std::string id;
	std::string seq;
	std::vector<size_t> polyApos;
	size_t len = 1;
};

void readKnownPolyA (const fs::path & f, std::vector<KnownPolyA> & knownPolyAvec);

void buildSeqStruct (std::vector<SeqStruct> & seqStt, const std::vector<KnownPolyA> & knownPolyAvec);

#endif /* __READKNOWNPOLYA_HPP__ */


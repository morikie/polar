#ifndef __POLARUTILITY_HPP__
#define __POLARUTILITY_HPP__

#include <string>
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <seqan/seq_io.h>
#include "refGeneParser.hpp" 

namespace polar {
namespace utility {

namespace fs = boost::filesystem;

size_t motifToIndex(std::string & motif);
char complement(const char c);
size_t getFastaIndex(const std::string & chr);
std::unordered_map<std::string, size_t> getTxRefSeqAccessions(const fs::path & f);
size_t getUtrLength(const RefGeneProperties & txProp);
std::string getUtrSequence(const RefGeneProperties & txProp, seqan::FaiIndex & fai);
size_t mapGenomePosToTxPos(const RefGeneProperties & txProp, size_t genomePos);
size_t getTxLength(const RefGeneProperties & txProp);
std::string getSeqAfterUtr(const RefGeneProperties & txProp, seqan::FaiIndex & faiIndex, size_t offset); 

} /* namespace utility */
} /* namespace polar */

#endif /* __POLARUTILITY_HPP__ */


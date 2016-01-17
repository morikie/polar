#ifndef __POLARUTILITY_HPP__
#define __POLARUTILITY_HPP__

#include <string>
#include <unordered_map>
#include <boost/filesystem.hpp>

namespace polar {
namespace utility {

namespace fs = boost::filesystem;

size_t motifToIndex(std::string & motif);
char complement(const char c);
size_t getFastaIndex(const std::string & chr);
std::unordered_map<std::string, size_t> getTxRefSeqAccessions(const fs::path & f);

} /* namespace utility */
} /* namespace polar */

#endif /* __POLARUTILITY_HPP__ */


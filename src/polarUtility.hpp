#ifndef __POLARUTILITY_HPP__
#define __POLARUTILITY_HPP__

#include<string>

namespace polar {
namespace utility {

size_t motifToIndex(std::string & motif);
char complement(const char c);
size_t getFastaIndex(const std::string & chr);

} /* namespace utility */
} /* namespace polar */

#endif /* __POLARUTILITY_HPP__ */


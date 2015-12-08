#include <string>
#include <stdexcept>
#include <boost/spirit/include/qi.hpp>
#include "polarUtility.hpp"

namespace polar {
namespace utility {

namespace qi = boost::spirit::qi;

size_t getFastaIndex(std::string & chr) {
	size_t idx = UINT_MAX;
	if (chr == "chrX") {
		idx = 22;
	} else if (chr == "chrY") {
		idx = 23;
	} else {
		if (! qi::parse(chr.begin(), chr.end(), qi::omit[*qi::alpha] >> qi::uint_, idx)) {
			throw std::invalid_argument("error parsing chromosome value from vcf");
	}
		//adjusting idx to 0-starting map
		idx--;
	}
	return idx; 
}

} /* namespace utility */
} /* namespace polar */


#include <climits>
#include <string>
#include <stdexcept>
#include <boost/spirit/include/qi.hpp>
#include "polarUtility.hpp"


namespace polar {
namespace utility {

namespace qi = boost::spirit::qi;


size_t motifToIndex(std::string & motif) {
	if (motif == "aataaa") return 0;
	else if (motif == "attaaa") return 1;
	else if (motif == "tataaa") return 2;
	else if (motif == "agtaaa") return 3;
	else if (motif == "aagaaa") return 4;
	else if (motif == "aatata") return 5;
	else if (motif == "aataca") return 6;
	else if (motif == "cataaa") return 7;
	else if (motif == "gataaa") return 8;
	else if (motif == "aatgaa") return 9;
	else if (motif == "actaaa") return 10;
	else if (motif == "aataga") return 11;
	else return 12;
}


char complement(const char c) {
	switch(c) {
	case 'a':
		return 't';
		break;
	case 'c':
		return 'g';
		break;
	case 'g':
		return 'c';
		break;
	case 't':
		return 'a';
		break;
	default:
		return 'n';
		break;
	}
}


size_t getFastaIndex(const std::string & chr) {
	size_t idx = UINT_MAX;
	if (chr == "chrX") {
		idx = 22;
	} else if (chr == "chrY") {
		idx = 23;
	} else {
		if (! qi::parse(chr.begin(), chr.end(), qi::omit[*qi::alpha] >> qi::uint_, idx)) {
			throw std::invalid_argument("error parsing chromosome for FAI");
	}
	//adjusting idx to 0-starting map
	idx--;
	}
	return idx; 
}

} /* namespace utility */
} /* namespace polar */


#include <climits>
#include <fstream>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/filesystem.hpp>
#include "polarUtility.hpp"


namespace polar {
namespace utility {

namespace qi = boost::spirit::qi;
namespace fs = boost::filesystem;

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
	if (chr == "chrX" || chr == "X" || chr == "x") {
		idx = 22;
	} else if (chr == "chrY" || chr == "Y" || chr == "y") {
		idx = 23;
	} else {
		if (! qi::parse(chr.begin(), chr.end(), qi::omit[*qi::alpha] >> qi::uint_, idx)) {
			return UINT_MAX; 
		}
		//adjusting idx to 0-starting map
		idx--;
	}
	return idx; 
}

std::unordered_map<std::string, size_t> getTxRefSeqAccessions(const fs::path & f) {
	typedef std::string accession;
	typedef size_t patch;
	std::unordered_map<accession, patch> txRefSeq;
	std::fstream in(f.string());
	std::string line;
	qi::rule<std::string::iterator, std::pair<accession, patch>() > r = *~qi::char_('\t') >> '\t' >> qi::uint_;
	while (std::getline(in, line)) {
		std::pair<accession, patch> temp;
		qi::parse(line.begin(), line.end(), r, temp);
		if (txRefSeq.find(temp.first) == txRefSeq.end()) txRefSeq[temp.first] = temp.second;
		else std::cerr << "Different patch numbers used for the same transcript" << std::endl;
	}
	return txRefSeq;
}

} /* namespace utility */
} /* namespace polar */


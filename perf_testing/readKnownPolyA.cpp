#include <algorithm>
#include <fstream>
#include <boost/filesystem/path.hpp>
#include <boost/spirit/include/qi.hpp>
#include "../src/hgvsParser.hpp"
#include "../src/seqStruct.hpp"
#include "readKnownPolyA.hpp"

namespace fs = boost::filesystem;
namespace qi = boost::spirit::qi;

void readKnownPolyA (const fs::path & f, std::vector<KnownPolyA> & knownPolyAvec) {
	std::ifstream in(f.string());
	std::string line;
	size_t lineCount = 0;
	KnownPolyA temp;

	while (std::getline(in, line)) {
		switch ( lineCount ) {
		case 0:
		{
			temp.id = std::string(line.begin() + 1, line.end());
			lineCount++;
			break;
		}
		case 1:
		{
			temp.polyApos.clear();
			qi::rule<std::string::iterator, std::vector<size_t>()> r = 
				qi::lit("Poly(A)_pos:") >> +(qi::uint_ >> qi::lit(','));
			qi::parse(line.begin(), line.end(), r, temp.polyApos);
			lineCount++;
			break;
		}
		case 2:
		{	
			std::transform(line.begin(), line.end(), line.begin(), ::tolower);
			temp.seq = line;
			knownPolyAvec.push_back(temp);
			lineCount = 0;
			break;
		}
		}
	}
}

void buildSeqStruct (std::vector<SeqStruct> & seqStt, const std::vector<KnownPolyA> & knownPolyAvec) {
	auto it = knownPolyAvec.begin();
	for (; it != knownPolyAvec.end(); it++) {
		SeqStruct seqSt = {
			it->seq,
			it->len,
			it->seq.size(),
			boost::none,
			boost::none,
			boost::none,
			boost::none,
			it->id
		};
		seqStt.push_back(seqSt);
	}
}



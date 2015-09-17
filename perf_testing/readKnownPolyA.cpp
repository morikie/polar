#include <fstream>
#include <boost/filesystem/path.hpp>
#include <boost/spirit/include/qi.hpp>
#include "../src/seqStruct.hpp"
#include "readKnownPolyA.hpp"

namespace fs = boost::filesystem;
namespace qi = boost::spirit::qi;

void readKnownPolyA (const fs::path & f, std::vector<KnownPolyA> & knownPolyAvec) {
	std::ifstream in(f.string());
	std::string line;
	unsigned int lineCount = 0;
	KnownPolyA temp;

	while (std::getline(in, line)) {
		switch ( lineCount ) {
		case 0:
		{
			temp.id = line;
			lineCount++;
			break;
		}
		case 1:
		{
			temp.polyApos.clear();
			qi::rule<std::string::iterator, std::vector<unsigned int>()> r = 
				qi::lit("Poly(A)_pos:") >> +(qi::uint_ >> qi::lit(','));
			qi::parse(line.begin(), line.end(), r, temp.polyApos);
			lineCount++;
			break;
		}
		case 2:
		{
			temp.seq = line;
			knownPolyAvec.push_back(temp);
			lineCount = 0;
			break;
		}
		}
	}
}

void buildSeqStruct (std::vector<SeqStruct> & seqStt, const std::vector<KnownPolyA> & knownPolyAvec) {
	auto itBegin = knownPolyAvec.begin();
	for (; itBegin != knownPolyAvec.end(); itBegin++) {
		SeqStruct transMut = {
			txSequences.getValueByKey(vcfTx.txName),
			utr3Start,
			txLength,
			boost::optional<const HgvsParser>(HgvsParser(vcfTx.hgvsString)),
			boost::optional<const std::string &>(chrom),
			boost::optional<const size_t &>(genePos),
			boost::optional<const std::string &>(strand),
			boost::optional<const std::string &>(vcfTx.txName)
			};
		seqStt.push_back( 
	}
}



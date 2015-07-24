#ifndef __KNOWNGENEPARSER_HPP__
#define __KNOWNGENEPARSER_HPP__

#include <string>
#include <unordered_map>
#include <boost/filesystem/path.hpp>


namespace fs = boost::filesystem;


struct TxProperties {

	friend bool operator== (const TxProperties & left, const TxProperties & right) {
		if (left.chr != right.chr) return false;
		if (left.strand != right.strand) return false;
		if (left.txStart != right.txStart) return false;
		if (left.txEnd != right.txEnd) return false;
		if (left.cdsStart != right.cdsStart) return false;
		if (left.cdsEnd != right.cdsEnd) return false;
		if (left.exonStarts != right.exonStarts) return false;
		if (left.exonEnds != right.exonEnds) return false;
		return true;
	}

	std::string chr;
	std::string strand;
	unsigned int txStart;
	unsigned int txEnd;
	unsigned int cdsStart;
	unsigned int cdsEnd;
	std::vector<unsigned int> exonStarts;
	std::vector<unsigned int> exonEnds;
};


class KnownGeneParser {
public:
	typedef std::unordered_map<std::string, TxProperties> txPropMap;

private:
	fs::path file;
	txPropMap data;
	
public:
	KnownGeneParser(const fs::path & f);
	~KnownGeneParser();

	txPropMap & getData();

private:
	void parse();

};

#endif /* __KNOWNGENEPARSER_HPP__ */


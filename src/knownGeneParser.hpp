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
	typedef std::string transcriptName;
	typedef std::unordered_map<transcriptName, TxProperties> txPropMap;

private:
	/* file path */
	fs::path file;
	
	/* empty TxProperty struct for return values */
	TxProperties emptyTxProperties = TxProperties { 
		std::string(),
		std::string(),
		0u,
		0u,
		0u,
		0u,
		std::vector<unsigned int> {},
		std::vector<unsigned int> {}
	};

	/* parsed data */
	txPropMap data;
	
public:
	KnownGeneParser(const fs::path & f);
	~KnownGeneParser();

	const TxProperties & getValueByKey(const std::string & k) const;

private:
	void parse();

};

#endif /* __KNOWNGENEPARSER_HPP__ */


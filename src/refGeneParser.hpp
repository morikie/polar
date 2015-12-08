#ifndef __REFGENEPARSER_HPP__
#define __REFGENEPARSER_HPP__

#include <string>
#include <unordered_map>
#include <boost/filesystem/path.hpp>


namespace fs = boost::filesystem;


struct RefGeneProperties {

	friend bool operator== (const RefGeneProperties & left, const RefGeneProperties & right) {
		if (left.chr != right.chr) return false;
		if (left.strand != right.strand) return false;
		if (left.txStart != right.txStart) return false;
		if (left.txEnd != right.txEnd) return false;
		if (left.cdsStart != right.cdsStart) return false;
		if (left.cdsEnd != right.cdsEnd) return false;
		if (left.exonStarts != right.exonStarts) return false;
		if (left.exonEnds != right.exonEnds) return false;
		if (left.score != right.score) return false;
		if (left.name2 != right.name2) return false;
		if (left.cdsStartStat != right.cdsStartStat) return false;
		if (left.cdsEndStat != right.cdsEndStat) return false;
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
	int score;
	std::string name2;
	std::string cdsStartStat;
	std::string cdsEndStat;

};


class RefGeneParser {
public:
	typedef std::string transcriptName;
	typedef std::unordered_map<transcriptName, RefGeneProperties> txPropMap;

	/* empty RefGeneProperies struct for return values */
	RefGeneProperties emptyRefGeneProperties = RefGeneProperties { 
		std::string(),
		std::string(),
		0u,
		0u,
		0u,
		0u,
		std::vector<unsigned int> {},
		std::vector<unsigned int> {},
		0,
		std::string(),
		std::string(),
		std::string()
	};

private:
	/* file path */
	fs::path file;
	/* parsed data */
	txPropMap data;
	
public:
	RefGeneParser(const fs::path & f);
	~RefGeneParser();

	const RefGeneProperties & getValueByKey(const std::string & k) const;
	std::vector<std::string> getKeys() const;

private:
	void parse();

};

#endif /* __REFGENEPARSER_HPP__ */


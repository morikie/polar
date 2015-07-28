#ifndef __JANNOVARVCFPARSER_HPP__
#define __JANNOVARVCFPARSER_HPP__

#include <string>
#include <map>
#include <boost/filesystem/path.hpp>


namespace fs = boost::filesystem;


/**
 * Struct containing desired values from the vcf file, i.e. variant type, transcript name
 * and HGVS string.
 */
struct vcfTranscripts {

	friend bool operator== (const vcfTranscripts & left, const vcfTranscripts & right) {
		if (left.jvVariantType != right.jvVariantType) return false;
		if (left.txName != right.txName) return false;
		if (left.hgvsString != right.hgvsString) return false;
		return true;
	}

	std::string jvVariantType;
	std::string txName;
	std::string hgvsString;
};


/**
 * Class to parse VCF files and fill vcfTranscripts.
 */
class JannovarVcfParser {
public:
	typedef std::string chromosome;
	typedef unsigned int position;
	typedef std::pair<chromosome, position> Key;
	typedef std::vector<vcfTranscripts> Value;
	typedef std::map<Key, Value> vcfTranscriptsMap;

private:
	fs::path file;
	vcfTranscriptsMap data;
	
public:
	JannovarVcfParser(const fs::path & f);
	~JannovarVcfParser();

	vcfTranscriptsMap & getData();

private:
	void parse();

};

#endif /* __JANNOVARVCFPARSER_HPP__ */


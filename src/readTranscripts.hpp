#ifndef __READTRANSCRIPTS_HPP__
#define __READTRANSCRIPTS_HPP__

#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <string>
#include <unordered_map>
#include <boost/filesystem/path.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>


namespace qi = boost::spirit::qi;
namespace fs = boost::filesystem;


/**
 * Class used to parse and store values from UCSC's knownGeneMrna files.
 */
class ReadTranscripts {	
public:
	typedef std::unordered_map<std::string, std::string> tMap;

private:
	/* file path */
	const fs::path file;
	
	/* empty string for return values */
	const std::string emptyString = std::string();
	
	/* the parsed transcripts */
	tMap transcripts;

public:
	ReadTranscripts(const fs::path & f);
	~ReadTranscripts();
	
	const std::string & getValueByKey(const std::string & k) const;

private:
	void parse();

};

#endif /* __READTRANSCRIPTS_HPP__ */


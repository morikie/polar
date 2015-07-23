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


class ReadTranscripts {	

public:
	typedef std::unordered_map<std::string, std::string> tMap;

private:
	const fs::path file;
	tMap transcripts;

public:
	ReadTranscripts(const fs::path & f);
	~ReadTranscripts();
	
	tMap & getData();

private:
	void parse();

};

#endif /* __READTRANSCRIPTS_HPP__ */

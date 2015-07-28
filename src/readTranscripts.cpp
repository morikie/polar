#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <fstream>
#include <string>
#include <unordered_map>
#include <boost/spirit/home/support/multi_pass.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include "readTranscripts.hpp"


namespace qi = boost::spirit::qi;
namespace spirit = boost::spirit;


/**
 * Grammar for parsing UCSC's knownGeneMrna files.
 */
template <typename Iterator>
struct knownGeneMrnaGrammar : 
	qi::grammar<Iterator, std::unordered_map<std::string, std::string>()> {
	
	/* Grammar */
	knownGeneMrnaGrammar() : knownGeneMrnaGrammar::base_type(query) {
			query   =  pair >> *(qi::lit('\n') >> pair);
			pair    =  seqName >> '\t' >> seq;
			seqName =  *~qi::char_('\t');
			seq     = +qi::char_("a-zA-Z_0-9");
	}

	qi::rule<Iterator, std::unordered_map<std::string, std::string>()> query;
	qi::rule<Iterator, std::pair<std::string, std::string>()> pair;
	qi::rule<Iterator, std::string()> seqName, seq;
};


/**
 * Constructor.
 *
 * @param[in]	f	File to be parsed.
 */
ReadTranscripts::ReadTranscripts(const fs::path & f):
	file(f)
{
	this->parse();	
}


/**
 * Destructor.
 */
ReadTranscripts::~ReadTranscripts() {}


/**
 * Return a reference to the parsed data.
 */
ReadTranscripts::tMap & ReadTranscripts::getData() {
	return this->transcripts;
}


/**
 * Parses the file passed to the constructor.
 */
void ReadTranscripts::parse() {
	std::ifstream in((this->file).string());

	typedef std::istreambuf_iterator<char> base_iterator_type;
	typedef spirit::multi_pass<base_iterator_type> forward_iterator;
	forward_iterator first = spirit::make_default_multi_pass(base_iterator_type(in));
	forward_iterator last = spirit::make_default_multi_pass(base_iterator_type());

	knownGeneMrnaGrammar<forward_iterator > grammar;

	
	bool result = qi::parse(first, last, grammar, this->transcripts);

}

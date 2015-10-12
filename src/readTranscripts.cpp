//#define BOOST_SPIRIT_DEBUG
#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <fstream>
#include <string>
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
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
	
	
	//BOOST_SPIRIT_DEBUG_NODES((query)(pair)(seq))
	}
private:
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
 * Returns the transcript sequence.
 *
 * @param[in]	k	Name of the transcript.
 */
const std::string & ReadTranscripts::getValueByKey(const std::string & k) const{
	auto it = this->transcripts.find(k);
	if (it != this->transcripts.end()) {
		return it->second;
	} else {
		return this->emptyString;
	}
}


/**
 * Parses the file passed to the constructor.
 */
void ReadTranscripts::parse() {
	if (! boost::filesystem::exists(this->file)) {
		std::cerr << "Could not find " << this->file.string() << std::endl;
		return;
	}
	boost::iostreams::mapped_file_source in(this->file.string());

	typedef char const* base_iterator_type;
	base_iterator_type first(in.data());
	base_iterator_type last(in.end());
	//typedef spirit::multi_pass<base_iterator_type> forward_iterator;
	//forward_iterator first = spirit::make_default_multi_pass(base_iterator_type(in));
	//forward_iterator last = spirit::make_default_multi_pass(base_iterator_type());

	knownGeneMrnaGrammar<base_iterator_type> grammar;

	
	bool result = qi::parse(first, last, grammar, this->transcripts);

}


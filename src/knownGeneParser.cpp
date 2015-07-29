#define BOOST_SPIRIT_DEBUG
#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <fstream>
#include <string>
#include <unordered_map>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/home/support/multi_pass.hpp>
#include "knownGeneParser.hpp"


namespace qi = boost::spirit::qi;
namespace spirit = boost::spirit;


BOOST_FUSION_ADAPT_STRUCT (
        TxProperties,
        (std::string, chr)
        (std::string, strand)
        (unsigned int, txStart)
        (unsigned int, txEnd)
        (unsigned int, cdsStart)
        (unsigned int, cdsEnd)
        (std::vector<unsigned int>, exonStarts)
        (std::vector<unsigned int>, exonEnds)
)


template <typename Iterator>
struct knownGeneGrammar : 
	qi::grammar<Iterator, std::unordered_map<std::string, TxProperties>() > {

	/* Grammar */
	knownGeneGrammar() : knownGeneGrammar::base_type(query) {
		query   %= pair >> *(qi::lit('\n') >> pair);
		pair    %= seqName >> txProps;
		seqName %= *~qi::char_('\t') >> '\t';
		intValue %= qi::uint_ >> '\t';
		intVector %= *(qi::uint_ >> ',') >> '\t';
		txProps	%= seqName
				>> seqName
				>> intValue
				>> intValue
				>> intValue
				>> intValue
				>> qi::omit[*~qi::char_('\t')] >> '\t'
				>> intVector
				>> intVector
				>> qi::omit[*~qi::char_('\n')]
				;
	
		//BOOST_SPIRIT_DEBUG_NODES( (query)(pair)(seqName)(txProps) )
	}
	private:
		qi::rule<Iterator, std::unordered_map<std::string, TxProperties>()> query;
		qi::rule<Iterator, std::pair<std::string, TxProperties>() > pair;
		qi::rule<Iterator, std::string()> seqName;
		qi::rule<Iterator, unsigned int()> intValue;
		qi::rule<Iterator, std::vector<unsigned int>() > intVector;
		qi::rule<Iterator, TxProperties()> txProps;
};


KnownGeneParser::KnownGeneParser(const fs::path & f) :
	file(f)
{
	this->parse();
}


KnownGeneParser::~KnownGeneParser() {}


const TxProperties & KnownGeneParser::getValueByKey(const std::string & k) const {
	auto it = this->data.find(k);
	if (it != this->data.end()) {
		return it->second;
	} else {
		return this->emptyTxProperties;
	}
}


void KnownGeneParser::parse() {
	std::ifstream in((this->file).string());

	typedef std::istreambuf_iterator<char> base_iterator_type;
	typedef spirit::multi_pass<base_iterator_type> forward_iterator;
	forward_iterator first = spirit::make_default_multi_pass(base_iterator_type(in));
	forward_iterator last = spirit::make_default_multi_pass(base_iterator_type());
	
	knownGeneGrammar<forward_iterator> grammar;

	bool result = qi::parse(first, last, grammar, this->data);

}


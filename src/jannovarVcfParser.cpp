#define BOOST_SPIRIT_DEBUG
#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <fstream>
#include <string>
#include <unordered_map>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/home/support/multi_pass.hpp>
#include "jannovarVcfParser.hpp"


namespace qi = boost::spirit::qi;
namespace spirit = boost::spirit;


BOOST_FUSION_ADAPT_STRUCT (
        vcfTranscripts,
        (std::string, jvVariantType)
        (std::string, txName)
        (std::string, hgvsString)
)


template <typename Iterator>
struct jannovarVcfGrammar : 
	qi::grammar<Iterator, JannovarVcfParser::vcfTranscriptsMap()> {

	/* Grammar */
	jannovarVcfGrammar() : jannovarVcfGrammar::base_type(query) {
		query   %= *headerLines >>
				pair >> *(qi::lit('\n') >> pair);
		
		headerLines = qi::omit[('#' >> *~qi::char_('\n') >> '\n')];
		pair    %= keyPair >> omittedCol
				>> omittedCol
				>> omittedCol
				>> omittedCol
				>> omittedCol
				>> jannovarStringValueVector;
		keyPair %= chromValue >> intValue;
		jannovarStringValueVector %= (*(jannovarStringValue - qi::char_(";")) >> ';' >> qi::omit[*~qi::char_('\n')]) | noAnnotation; 
		jannovarStringValue %=	qi::omit[stringValue] >> '|'	
					>> stringValue >> '|' 
					>> qi::omit[stringValue] >> '|'
					>> qi::omit[stringValue] >> '|' 
					>> qi::omit[stringValue] >> '|' 
					>> qi::omit[stringValue] >> '|'
					>> stringValue >> '|'
					>> qi::omit[stringValue] >> '|'
					>> qi::omit[stringValue] >> '|'
					>> stringValue >> '|'
					>> qi::omit[stringValue] >> '|'
					>> qi::omit[stringValue] >> '|'
					>> qi::omit[stringValue] >> '|'
					>> qi::omit[stringValue] >> '|'
					>> qi::omit[stringValue] >> '|'
					;
		noAnnotation = qi::attr("") >> qi::attr("") >> qi::attr("") >> qi::omit[*~qi::char_('\n')];
		omittedCol = qi::omit[*~qi::char_('\t')] >> '\t';
		stringValue %= *~qi::char_("|\t");
		intValue %= qi::uint_ >> '\t';
		chromValue %= *~qi::char_('\t') >> '\t';
	
		BOOST_SPIRIT_DEBUG_NODES( (query)(pair)(keyPair)(headerLines)(jannovarStringValue)(jannovarStringValueVector)(stringValue) )
	}
	private:
		qi::rule<Iterator, JannovarVcfParser::vcfTranscriptsMap()> query;
		qi::rule<Iterator, JannovarVcfParser::Key()> keyPair;
		qi::rule<Iterator, JannovarVcfParser::Value()> jannovarStringValueVector;
		qi::rule<Iterator, std::pair<JannovarVcfParser::Key, JannovarVcfParser::Value>() > pair;
		qi::rule<Iterator, JannovarVcfParser::chromosome()> stringValue;
		qi::rule<Iterator, JannovarVcfParser::position()> intValue;
		qi::rule<Iterator, std::string()> chromValue;
		qi::rule<Iterator, vcfTranscripts()> jannovarStringValue;
		qi::rule<Iterator, vcfTranscripts()> noAnnotation;
		qi::rule<Iterator, void()> omittedCol;
		qi::rule<Iterator, void()> headerLines;
};


JannovarVcfParser::JannovarVcfParser(const fs::path & f) :
	file(f)
{
	this->parse();
}


JannovarVcfParser::~JannovarVcfParser() {}


JannovarVcfParser::vcfTranscriptsMap & JannovarVcfParser::getData() {
	return this->data;
}


void JannovarVcfParser::parse() {
	std::ifstream in((this->file).string());

	typedef std::istreambuf_iterator<char> base_iterator_type;
	typedef spirit::multi_pass<base_iterator_type> forward_iterator;
	forward_iterator first = spirit::make_default_multi_pass(base_iterator_type(in));
	forward_iterator last = spirit::make_default_multi_pass(base_iterator_type());
	
	jannovarVcfGrammar<forward_iterator> grammar;

	bool result = qi::parse(first, last, grammar, this->data);

}


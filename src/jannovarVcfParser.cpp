//#define BOOST_SPIRIT_DEBUG 
#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <fstream>
#include <string>
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/home/support/multi_pass.hpp>
#include "jannovarVcfParser.hpp"


namespace qi = boost::spirit::qi;
namespace spirit = boost::spirit;


/**
 *  boost::fusion "struct" of vcfTranscripts struct defined in the header.
 */
BOOST_FUSION_ADAPT_STRUCT (
        vcfTranscripts,
        (std::string, jvVariantType)
        (std::string, txName)
        (std::string, hgvsString)
)


/**
 * Grammar struct used by the parse function.
 */
template <typename Iterator>
struct jannovarVcfGrammar : 
	qi::grammar<Iterator, JannovarVcfParser::vcfTranscriptsMap()> {

	jannovarVcfGrammar() : jannovarVcfGrammar::base_type(query) {
		
		
		/* Grammar elements */	
		stringColumnValue %= *~qi::char_('\t') >> '\t';			//parser for tab seperated columns
		headerLines = ('#' >> *~qi::char_('\n') >> '\n');		//ignore header lines at start of file
		intColumnValue %= qi::uint_ >> '\t';				
		noAnnotation = qi::attr("") >> qi::attr("") >> qi::attr("");	//dummy values if jannovar string cannot be parsed (fills the fusion struct)
		stringValue %= *~qi::char_("|\t");				//parse chars until you encounter '|' or '\t'	
		
		/* Actual grammar structure. query being the entry point. */
		query	%= qi::omit[*headerLines]
			>> pair >> *(qi::eol >> pair);	//pair as in std::map<key, kalue>

		pair    %= keyPair  			//CHROM and POS colums (used as key for the std::map)
			>> qi::omit[stringColumnValue]	//ID column
			>> qi::omit[stringColumnValue]	//REF colum
			>> qi::omit[stringColumnValue] 	//ALT colum
			>> qi::omit[stringColumnValue] 	//QUAL colum
			>> qi::omit[stringColumnValue]	//FILTER colum
			>> qi::omit[*(qi::char_ - qi::lit("ANN="))] >> jannovarStringValueVector	//INFO column (contains the jannovar annotation string)
			>> qi::omit[*(qi::char_ - qi::eol)];	//skip anything that follows after the INFO column until EOL
		
		keyPair %= stringColumnValue >> intColumnValue;
		
		//parse any number of jannovar strings or stop parsing (noAnnotation) the line if '\t' is encountered (see stringValue) 
		jannovarStringValueVector %= (*(jannovarStringValue - qi::char_(";")) >> ';') | noAnnotation;
		//jannovar string consists of 15 columns. Three values are extracted (variant type, transcript and hgvs string) 
		jannovarStringValue %=	
			qi::omit[stringValue] >> '|'
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
	
	//BOOST_SPIRIT_DEBUG_NODES( (query)(pair)(keyPair)(headerLines)(jannovarStringValue)(jannovarStringValueVector)(stringValue) )
	}

private:
	qi::rule<Iterator, void()> headerLines;
	qi::rule<Iterator, JannovarVcfParser::position()> intColumnValue;
	qi::rule<Iterator, vcfTranscripts()> jannovarStringValue;
	qi::rule<Iterator, JannovarVcfParser::Value()> jannovarStringValueVector;
	qi::rule<Iterator, JannovarVcfParser::Key()> keyPair;
	qi::rule<Iterator, vcfTranscripts()> noAnnotation;
	qi::rule<Iterator, std::pair<JannovarVcfParser::Key, JannovarVcfParser::Value>() > pair;
	qi::rule<Iterator, JannovarVcfParser::vcfTranscriptsMap()> query;
	qi::rule<Iterator, JannovarVcfParser::chromosome()> stringColumnValue;
	qi::rule<Iterator, std::string()> stringValue;
};


/**
 * Constructor.
 *
 * @param[in]	f 	File to be parsed.
 */
JannovarVcfParser::JannovarVcfParser(const fs::path & f) :
	file(f)
{
	this->parse();
}


/**
 * Destructor. 
 */
JannovarVcfParser::~JannovarVcfParser() {}


/**
 * Returns a reference to the parsed data. 
 */
const JannovarVcfParser::vcfTranscriptsMap & JannovarVcfParser::getData() const {
	return this->data;
}


/**
 * Parses the file passed to the constructor. 
 */
void JannovarVcfParser::parse() {
	if (! boost::filesystem::exists(this->file)) {
		std::cerr << "Could not find " << this->file.string() << std::endl;
		return;
	}
	boost::iostreams::mapped_file_source  in(this->file.string());
	
	typedef char const* base_iterator_type;
	base_iterator_type first(in.data());
	base_iterator_type last(in.end());
	//typedef spirit::multi_pass<base_iterator_type> forward_iterator;
	//forward_iterator first = spirit::make_default_multi_pass(base_iterator_type(in));
	//forward_iterator last = spirit::make_default_multi_pass(base_iterator_type());
	
	jannovarVcfGrammar<base_iterator_type> grammar;

	bool result = qi::parse(first, last, grammar, this->data);
}


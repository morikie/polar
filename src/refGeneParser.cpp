//#define BOOST_SPIRIT_DEBUG
#define BOOST_SPIRIT_USE_PHOENIX_V3
#define FUSION_MAX_VECTOR_SIZE 15

#include <fstream>
#include <string>
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/home/support/multi_pass.hpp>
#include "refGeneParser.hpp"


namespace qi = boost::spirit::qi;
namespace spirit = boost::spirit;


BOOST_FUSION_ADAPT_STRUCT (
        RefGeneProperties,
        (std::string, chr)
        (std::string, strand)
        (unsigned int, txStart)
        (unsigned int, txEnd)
        (unsigned int, cdsStart)
        (unsigned int, cdsEnd)
        (std::vector<unsigned int>, exonStarts)
        (std::vector<unsigned int>, exonEnds)
	(int, score)
	(std::string, name2)
	(std::string, cdsStartStat)
	(std::string, cdsEndStat)
)


template <typename Iterator>
struct refGeneGrammar : 
	qi::grammar<Iterator, std::unordered_map<std::string, RefGeneProperties>() > {

	/* Grammar */
	refGeneGrammar() : refGeneGrammar::base_type(query) {
		query   %= pair >> *(qi::lit('\n') >> pair);
		pair    %= qi::omit[uIntValue] >> stringValue >> txProps;
		stringValue %= *~qi::char_('\t') >> '\t';
		intValue %= qi::int_ >> '\t';
		uIntValue %= qi::uint_ >> '\t';
		uIntVector %= *(qi::uint_ >> ',') >> '\t';
		txProps	%= stringValue
				>> stringValue
				>> uIntValue
				>> uIntValue
				>> uIntValue
				>> uIntValue
				>> qi::omit[*~qi::char_('\t')] >> '\t'
				>> uIntVector
				>> uIntVector
				>> intValue
				>> stringValue
				>> stringValue
				>> stringValue
				>> qi::omit[*~qi::char_('\n')]
		;
	
	//BOOST_SPIRIT_DEBUG_NODES( (query)(pair)(stringValue)(uIntValue)(uIntVector)(txProps))
	}

	private:
		qi::rule<Iterator, std::unordered_map<std::string, RefGeneProperties>()> query;
		qi::rule<Iterator, std::pair<std::string, RefGeneProperties>() > pair;
		qi::rule<Iterator, std::string()> stringValue;
		qi::rule<Iterator, int()> intValue;
		qi::rule<Iterator, unsigned int()> uIntValue;
		qi::rule<Iterator, std::vector<unsigned int>() > uIntVector;
		qi::rule<Iterator, RefGeneProperties()> txProps;
};


RefGeneParser::RefGeneParser(const fs::path & f) :
	file(f)
{
	this->parse();
}


RefGeneParser::~RefGeneParser() {}


const RefGeneProperties & RefGeneParser::getValueByKey(const std::string & k) const {
	auto it = this->data.find(k);
	if (it != this->data.end()) {
		return it->second;
	} else {
		return this->emptyRefGeneProperties;
	}
}


void RefGeneParser::parse() {
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
	
	refGeneGrammar<base_iterator_type> grammar;

	bool result = qi::parse(first, last, grammar, this->data);

}


#define BOOST_SPIRIT_DEBUG

#include <boost/spirit/include/qi.hpp>
#include "hgvsParser.hpp"


namespace qi = boost::spirit::qi;


/**
 * Constructor.
 */
HgvsParser::HgvsParser(const std::string & hgvs) :
	utr3MutPos(-1),
	hgvsString(hgvs)
{
	this->parse();
}


/**
 * Destructor.
 */
HgvsParser::~HgvsParser() {}


/**
 * Parsing the HGVS string passed to the constructor.
 */
void HgvsParser::parse() {
	auto itBegin = this->hgvsString.begin();
	auto itEnd = this->hgvsString.end();

	qi::rule<std::string::iterator, unsigned int()> r = qi::omit[*~qi::char_('*')] >> '*' >> qi::uint_ >> qi::lit("+-") | qi::omit[*qi::char_];
	//BOOST_SPIRIT_DEBUG_NODES( (r) )
	qi::parse(itBegin, itEnd, r, this->utr3MutPos);

}


int getMutPosition() {
	return this->utr3MutPos;
}


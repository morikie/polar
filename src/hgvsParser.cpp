#include <boost/spirit/include/qi.hpp>
#include "hgvsParser.hpp"


namespace qi = boost::spirit::qi;


/**
 *
 */
HgvsParser::HgvsParser(const std::string & hgvs) :
	utr3MutPos(0u),
	hgvsString(hgvs)
{
	this->parse();
}


/**
 *
 */
HgvsParser::~HgvsParser() {}


/**
 *
 */
void HgvsParser::parse() {
	auto itBegin = this->hgvsString.begin();
	auto itEnd = this->hgvsString.end();

	qi::rule<std::string::iterator, unsigned int> r = qi::omit[*~qi::char_('*')] >> '*' >> qi::uint_ >> qi::omit[*qi::char_];
	qi::parse(itBegin, itEnd, r, this->utr3MutPos);
}


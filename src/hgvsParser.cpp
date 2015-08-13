#define BOOST_SPIRIT_DEBUG

#include <boost/phoenix.hpp>
#include <boost/spirit/home/classic/actor/assign_actor.hpp>
#include <boost/spirit/include/qi.hpp>
#include "hgvsParser.hpp"


using boost::phoenix::ref;

namespace qi = boost::spirit::qi;
namespace classic = boost::spirit::classic;


/**
 * Constructor.
 */
HgvsParser::HgvsParser(const std::string & hgvs) :
	hgvsString(hgvs),
	utr3MutPos(-1)
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
	
	qi::rule<std::string::iterator, void()> notExonic = ((qi::lit('+') | qi::lit('-') | qi::lit('_')) >> qi::omit[*qi::char_]); 
	qi::rule<std::string::iterator, void()> exonic = (qi::omit[*qi::char_]);
	qi::rule<std::string::iterator, void()> entry = qi::omit[*~qi::char_('*')] >> '*' >> qi::int_[boost::phoenix::ref(this->utr3MutPos) = qi::_1 - 1] 
		>> notExonic[classic::assign_a(this->intronic, true)] | exonic[classic::assign_a(this->intronic, false)];
	//BOOST_SPIRIT_DEBUG_NODES( (r) )
	qi::parse(itBegin, itEnd, entry, this->utr3MutPos);
}


/**
 * Return the position of the mutation mapped to index starting with zero (hence the -1).
 */
int HgvsParser::getMutPosition() const {
	return this->utr3MutPos;
}


/**
 *
 */
bool HgvsParser::isIntronic() const {
	return this->intronic;
}

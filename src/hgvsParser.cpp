#include <climits>
#include <boost/phoenix.hpp>
#include <boost/spirit/home/classic/actor/assign_actor.hpp>
#include <boost/spirit/include/qi.hpp>
#include "hgvsParser.hpp"


namespace qi = boost::spirit::qi;
namespace classic = boost::spirit::classic;


/**
 * Constructor.
 */
HgvsParser::HgvsParser(const std::string & hgvs) :
	hgvsString(hgvs),
	utr3MutPos(INT_MIN)
{
	this->parse();
}


/**
 * Destructor.
 */
HgvsParser::~HgvsParser() {}


/**
 * Parsing the HGVS string passed to the constructor.
 *
 * TODO: more sophisticated Parser, maybe define grammar struct. (detection limited to substitutions/deletions)
 */
void HgvsParser::parse() {
	auto itBegin = this->hgvsString.begin();
	auto itEnd = this->hgvsString.end();
	
	qi::rule<std::string::iterator, void()> notExonic = (qi::lit('+') | qi::lit('-') | qi::lit('_')) >> qi::omit[*qi::char_]; 
	qi::rule<std::string::iterator, void()> deletion = qi::lit('_') >> qi::lit('*') >> qi::omit[qi::int_] >> qi::lit("del") >> qi::omit[*qi::char_];
	qi::rule<std::string::iterator, void()> exonic = qi::omit[*qi::char_];
	qi::rule<std::string::iterator, void()> entry = qi::omit[*~qi::char_('*')] >> qi::lit('*') >> qi::int_[boost::phoenix::ref(this->utr3MutPos) = qi::_1 - 1]
		>> ( 
			deletion[classic::assign_a(this->intronic, false)] | 
			notExonic[classic::assign_a(this->intronic, true)] | 
			exonic[classic::assign_a(this->intronic, false)]
	);
	if (! qi::parse(itBegin, itEnd, entry, this->utr3MutPos)) {
		this->utr3MutPos = INT_MIN;
		this->intronic = true;
	}
}


/**
 * Return the position of the mutation mapped to index starting with zero.
 */
int HgvsParser::getMutPosition() const {
	return this->utr3MutPos;
}


/**
 * Returns if variant is intronic.
 */
bool HgvsParser::isIntronic() const {
	return this->intronic;
}


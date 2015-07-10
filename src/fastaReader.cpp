//fastaReader.cpp

#include <iomanip>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/spirit/include/classic_position_iterator.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/qi.hpp>
#include "fastaReader.hpp"


using namespace std;

namespace fs = boost::filesystem;
namespace qi = boost::spirit::qi;
namespace pt = boost::posix_time;

template <typename Iterator, typename Skipper>
struct FastaGrammar : qi::grammar<Iterator, FastaReader::fastaVector(), qi::locals<string>, Skipper> {
	qi::rule<Iterator> infoLineStart;
	qi::rule<Iterator> inputEnd;
	qi::rule<Iterator> lineEnd;
	qi::rule<Iterator, string(), Skipper> infoLine;
	qi::rule<Iterator, string(), Skipper> seqLine;
	qi::rule<Iterator, FastaReader::fastaVector(), qi::locals<string>, Skipper> fasta;


	FastaGrammar() : FastaGrammar::base_type(fasta, "fasta") {
		using boost::spirit::standard::char_;
		using boost::phoenix::bind;
		using qi::eoi;
		using qi::eol;
		using qi::eps;
		using qi::lexeme;
		using qi::lit;
		using qi::_1;
		using qi::_val;
		using namespace qi::labels;
		
		infoLineStart = char_('>');
		inputEnd = eoi;
		
		/* grammar elements */		
		infoLine = lexeme[*(char_ - eol)];
		seqLine = *(char_ - infoLineStart);

		/* the grammar */
		fasta = *(infoLineStart > infoLine[_a = _1] 
			> seqLine[bind(&FastaGrammar::addValue, _val, _a, _1)]
			)
			> inputEnd
		;
		
		/* define sub-grammar names for syntax error messages */
		infoLineStart.name(">");
		infoLine.name("sequence identifier");
		seqLine.name("sequence");
		
	}

	/**
	 * The method is being used by class FastaReader to add
	 * info/sequence pairs to the passed FastaReader::fastaVector argument.
	 *
	 * @param[out] fa - add element to this object
	 * @param[in] info - element key to add
	 * @param[in] seq - element value to add
	 */
	static void addValue(FastaReader::fastaVector & fa, const string & info, const string & seq) {
		fa.push_back(make_pair(info, seq));
	}
};


FastaReader::FastaReader(const fs::path & f) {
	this->file = f;	
	this->parse();
}


FastaReader::~FastaReader() {}


const fs::path & FastaReader::getFile() const {
	return this->file;
}


const FastaReader::fastaVector::const_iterator FastaReader::getBeginIterator() const {
	return this->fV.cbegin();
}


const FastaReader::fastaVector::const_iterator FastaReader::getEndIterator() const {
	return this->fV.cend();
}


void FastaReader::parse() {
	if ( this->file.empty() ) throw string("FastaReader: No file specific.");
	if ( ! fs::is_regular_file(this->file) ) throw (string("FastaReader: File not found: ") + this->file.string());

	typedef boost::spirit::istream_iterator iterator_type;
	typedef boost::spirit::classic::position_iterator2<iterator_type> pos_iterator_type;
	typedef FastaGrammar<pos_iterator_type, boost::spirit::ascii::space_type> fastaGr;
	
	
	std::cerr << "Measuring: Read-in." << std::endl;
	const pt::ptime startMeasurement = pt::microsec_clock::universal_time();

	fs::ifstream fin(this->file);
	if ( ! fin.is_open() ) {
		throw (string("FastaReader: Access denied for: ") + this->file.string());
	}

	const pt::ptime endMeasurement = pt::microsec_clock::universal_time();
	pt::time_duration duration(endMeasurement - startMeasurement);
	std::cerr << duration << std::endl;

	fin.unsetf(ios::skipws);

	/* input iterator as forward iterator, usable by spirit parser */
	iterator_type begin(fin);
	iterator_type end;

	/* wrap forward iterator with position iterator, to record the position */
	pos_iterator_type pos_begin(begin, end, this->file.string());
	pos_iterator_type pos_end;

	fastaGr fG;
	try {
		std::cerr << "Measuring: Parsing." << std::endl;
		const pt::ptime startMeasurement2 = pt::microsec_clock::universal_time();
		
		qi::phrase_parse(pos_begin, pos_end, fG, boost::spirit::ascii::space, this->fV);
			
		const pt::ptime endMeasurement2 = pt::microsec_clock::universal_time();
		pt::time_duration duration2(endMeasurement2 - startMeasurement2);
		std::cerr << duration2 <<  std::endl;
	} catch (std::string str) {
		cerr << "error message: " << str << endl;
	}	
}


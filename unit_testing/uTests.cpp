#define BOOST_TEST_MAIN polar 
#include <boost/test/unit_test.hpp>
#include "../src/polar.hpp"
#include "../src/fastaReader.hpp"

namespace fs = boost::filesystem;


BOOST_AUTO_TEST_CASE( readIn )
{
	fs::path f = "test.txt";
	FastaReader newReader(f);
}



//fastaRreader.hpp
#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;

class FastaReader {
public:
	FastaReader(const fs::path & f);
};


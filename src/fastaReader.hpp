//fastaRreader.hpp
#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;


class FastaReader {

public:
	typedef std::vector< std::pair<std::string, std::string> > fastaVector;

private:
	fastaVector fV;
	fs::path file;	

public:
	FastaReader(const fs::path & f);
	~FastaReader();

	const fs::path & getFile() const;
	const fastaVector::const_iterator getBeginIterator() const;
	const fastaVector::const_iterator getEndIterator() const;	

private:
	void parse();
	
};


#include <boost/filesystem/path.hpp>
#include <seqan/seq_io.h>
#include "buildIndexFile.hpp"

namespace fs = boost::filesystem;


void buildIndexFile(const fs::path & f) {
	seqan::FaiIndex faiIndex;
	if (!seqan::build(faiIndex, f.string().c_str())) {
		std::cerr << "error: could not build index file" << std::endl;
		return;
	}
	std::string indexPath = f.string() + ".fai";
	if (!seqan::save(faiIndex, indexPath.c_str() )) {
		std::cerr << "error: could not save index file" << std::endl;
		return;
	}
}


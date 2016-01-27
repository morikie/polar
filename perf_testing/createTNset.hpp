#ifndef __CREATETNSET_HPP__
#define __CREATETNSET_HPP__

#include <boost/filesystem.hpp>
#include <seqan/seq_io.h>

namespace fs = boost::filesystem;

bool createTNset(const fs::path & out, const seqan::FaiIndex & refGenomeIndex);

#endif /* __CREATETNSET_HPP__ */

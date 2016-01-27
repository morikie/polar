#ifndef __CREATETPSET_HPP__
#define __CREATETPSET_HPP__

#include <boost/filesystem.hpp>
#include <seqan/seq_io.h>
namespace fs = boost::filesystem;

bool createTPset(const fs::path & out, const seqan::FaiIndex & refGenomeIndex);

#endif /*__CREATETPSET_HPP__ */

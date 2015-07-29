#ifndef __HGVSPARSER_HPP__
#define __HGVSPARSER_HPP__

#include <string>


class HgvsParser {
public:
	unsigned int utr3MutPos;

private:
	std::string hgvsString;
	
	/* position of the mutation in the utr3 (after stop codon) */

public:
	HgvsParser(const std::string & hgvs);
	~HgvsParser();

private:
	void parse();

};

#endif /* __HGVSPARSER_HPP__ */


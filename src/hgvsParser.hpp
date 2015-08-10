#ifndef __HGVSPARSER_HPP__
#define __HGVSPARSER_HPP__

#include <string>


class HgvsParser {
private:
	std::string hgvsString;
	/* position of the mutation in the utr3 (after stop codon) */
	int utr3MutPos;
	bool intronic;

public:
	HgvsParser(const std::string & hgvs);
	~HgvsParser();

	int getMutPosition() const;
	bool isIntronic() const;

private:
	void parse();

};

#endif /* __HGVSPARSER_HPP__ */


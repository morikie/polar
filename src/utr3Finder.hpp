#ifndef __UTR3FINDER_HPP__
#define __UTR3FINDER_HPP__

#include <string>
#include <vector>
#include "seqStruct.hpp"

class Utr3Finder {
public:
	static const std::vector<std::string> hexamers;
	static const std::vector<std::string> rHexamers;

	static const size_t noHitPos = -1;

protected:
	const SeqStruct & seqStruct;

public:
	Utr3Finder(const SeqStruct & sSt);
	virtual ~Utr3Finder() = 0;

	virtual size_t getPolyaMotifPos() const  = 0;
	virtual bool isMutationInMotif() const = 0;
	virtual void writeInfo() const  = 0;

protected:
	virtual void findPolyaMotif() = 0;

};

#endif /* __UTR3FINDER_HPP__ */

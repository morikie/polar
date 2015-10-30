#ifndef __UTR3FINDER_HPP__
#define __UTR3FINDER_HPP__

#include <string>
#include <vector>
#include "seqStruct.hpp"

/* Interface class */
class Utr3Finder {
public:
	static const std::vector<std::string> hexamers;
	static const std::vector<std::string> rHexamers;

	static const size_t noHitPos = -1;

protected:
	const SeqStruct & seqStruct;
	std::vector<size_t> polyaPosVector = std::vector<size_t>(1, noHitPos);
	
public:
	Utr3Finder(const SeqStruct & sSt);
	virtual ~Utr3Finder() = 0;

	virtual std::vector<size_t> getPolyaMotifPos() const  = 0;
	virtual std::string getSequence() const = 0;
	virtual std::string getMotifSequence(const size_t & pos) const = 0;
	virtual bool isMutationInMotif() const = 0;
	virtual void writeInfo() const  = 0;

protected:
	virtual void findPolyaMotif() = 0;

};

#endif /* __UTR3FINDER_HPP__ */

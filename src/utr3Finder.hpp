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

	struct Utr3FinderResult {
		size_t pos;
		double truthValue;
		std::string strand;
	};

protected:
	const SeqStruct & seqStruct;
	std::vector<Utr3FinderResult> polyaPosVector;
	
public:
	Utr3Finder(const SeqStruct & sSt);
	virtual ~Utr3Finder() = 0;

	virtual std::vector<Utr3FinderResult> getPolyaMotifPos() const  = 0;
	virtual std::string getSequence() const = 0;
	virtual std::string getMotifSequence(const Utr3FinderResult & result) const = 0;
	virtual bool isMutationInMotif() const = 0;
	virtual void writeInfo() const  = 0;
	
protected:
	virtual void findPolyaMotif() = 0;

};

#endif /* __UTR3FINDER_HPP__ */

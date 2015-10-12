#ifndef __UTR3FINDERNAIVE_HPP__
#define __UTR3FINDERNAIVE_HPP__

#include "seqStruct.hpp"
#include "utr3Finder.hpp"


/**
 * Class to predict the location of the Poly(A) cleavage motif.
 */
class Utr3FinderNaive : public Utr3Finder {	
protected:
	size_t polyaMotifPos = noHitPos;
	
public:
	Utr3FinderNaive(const SeqStruct & sSt);
	virtual ~Utr3FinderNaive();

	virtual bool isMutationInMotif() const;
	virtual std::string getSequence() const;
	virtual std::string getMotifSequence() const;
	virtual size_t getPolyaMotifPos() const;
	virtual void writeInfo() const;

protected:
	virtual void findPolyaMotif();
	
};

#endif /* __UTR3FINDERNAIVE_HPP__ */


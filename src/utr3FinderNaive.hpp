#ifndef __UTR3FINDERNAIVE_HPP__
#define __UTR3FINDERNAIVE_HPP__

#include <string>
#include <vector>
#include "seqStruct.hpp"
#include "utr3Finder.hpp"


/**
 * Class to predict the location of the Poly(A) cleavage motif.
 */
class Utr3FinderNaive : public Utr3Finder {	
public:
	Utr3FinderNaive(const SeqStruct & sSt);
	virtual ~Utr3FinderNaive();

	virtual bool isMutationInMotif() const override;
	virtual std::string getSequence() const override;
	virtual std::vector<std::string> getMotifSequence() const override;
	virtual std::vector<size_t> getPolyaMotifPos() const override;
	virtual void writeInfo() const override;

protected:
	virtual void findPolyaMotif() override;
	
};

#endif /* __UTR3FINDERNAIVE_HPP__ */


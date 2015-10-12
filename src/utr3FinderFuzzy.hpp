#ifndef __UTR3FINDERFUZZY_HPP__
#define __UTR3FINDERFUZZY_HPP__

#include <string>
#include "seqStruct.hpp"
#include "utr3Finder.hpp"


class Utr3FinderFuzzy : public Utr3Finder {
protected:
	size_t polyaMotifPos = noHitPos;

public: 
	Utr3FinderFuzzy(const SeqStruct & sSt);
	virtual ~Utr3FinderFuzzy();

	bool isMutationInMotif() const;
	std::string getSequence() const;
	std::string getMotifSequence() const;
	size_t getPolyaMotifPos() const;
	void writeInfo() const;

protected:
	virtual void findPolyaMotif();

};      

#endif /* __UTR3FINDERFUZZY_HPP__ */

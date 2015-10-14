#ifndef __UTR3FINDERFUZZY_HPP__
#define __UTR3FINDERFUZZY_HPP__

#include <string>
#include "seqStruct.hpp"
#include "utr3Finder.hpp"


class Utr3FinderFuzzy : public Utr3Finder {
public: 
	Utr3FinderFuzzy(const SeqStruct & sSt);
	virtual ~Utr3FinderFuzzy();

	virtual bool isMutationInMotif() const override;
	virtual std::string getSequence() const override;
	virtual std::vector<std::string> getMotifSequence() const override;
	virtual std::vector<size_t> getPolyaMotifPos() const override;
	virtual void writeInfo() const override;

protected:
	virtual void findPolyaMotif() override;

};      

#endif /* __UTR3FINDERFUZZY_HPP__ */

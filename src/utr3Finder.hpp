#ifndef __UT3MUTATIONFINDER_HPP__
#define __UT3MUTATIONFINDER_HPP__

#include <boost/optional.hpp>
#include "seqStruct.hpp"


/**
 * Class to predict the location of the Poly(A) cleavage motif.
 */
class Utr3Finder {
public:
	static const std::vector<std::string> hexamers;
	static const std::vector<std::string> rHexamers;
	
	static const size_t noHitPos = -1;

private:
	const SeqStruct & txMut;
	size_t polyaMotifPos = noHitPos;
	
public:
	Utr3Finder(const SeqStruct & tM);
	~Utr3Finder();

	bool isMutationInMotif() const;
	std::string getSequence() const;
	std::string getMotifSequence() const;
	size_t getPolyaMotifPos() const;
	std::string writeInfo() const;

private:
	void findPolyaMotif();
	
};

#endif /* __UTR3MUTATIONFINDER_HPP__ */


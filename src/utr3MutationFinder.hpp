#ifndef __UT3MUTATIONFINDER_HPP__
#define __UT3MUTATIONFINDER_HPP__

#include "transcriptMutation.hpp"


/**
 * Class to predict the location of the Poly(A) cleavage motif.
 */
class Utr3MutationFinder {
public:
	static const std::vector<std::string> hexamers;
//	const TranscriptMutation mutation;

private:
	const TranscriptMutation & txMut;
	unsigned int polyaMotifPos;

public:
	Utr3MutationFinder(const TranscriptMutation & tM);
	~Utr3MutationFinder();

	bool isMutationInMotif() const;

private:
	void findPolyaMotif();
	
};

#endif /* __UTR3MUTATIONFINDER_HPP__ */


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
	static const std::vector<std::string> rHexamers;
	
	static const size_t noHitPos = -1;

	const TranscriptMutation & txMut;
private:
	size_t polyaMotifPos = 0;

public:
	Utr3MutationFinder(const TranscriptMutation & tM);
	~Utr3MutationFinder();

	bool isMutationInMotif() const;
	std::string getMotifSequence() const;
	size_t getPolyaMotifPos() const;
	std::string writeLocation() const;

private:
	void findPolyaMotif();
	
};

#endif /* __UTR3MUTATIONFINDER_HPP__ */


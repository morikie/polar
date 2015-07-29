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
	const std::string & sequence;

public:
	Utr3MutationFinder(const std::string & seq);
	~Utr3MutationFinder();

	void findConsensus();

};

#endif /* __UTR3MUTATIONFINDER_HPP__ */


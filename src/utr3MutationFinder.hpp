
#ifndef __UT3MUTATIONFINDER_HPP__
#define __UT3MUTATIONFINDER_HPP__

#include "transcriptMutation.hpp"


class Utr3MutationFinder {
public:
	static const std::set<std::string> hexamers;
//	const TranscriptMutation mutation;

private:
	const std::string & sequence;

public:
	Utr3MutationFinder(const std::string & seq);
	~Utr3MutationFinder();

	void findConsensus();

};

#endif /* __UTR3MUTATIONFINDER_HPP__ */


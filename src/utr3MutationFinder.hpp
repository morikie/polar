
class MutationFinder {
public:
	static const std::set<std::string> hexamers;

private:
	const std::string & sequence;

public:
	MutationFinder(const std::string & seq);
	~MutationFinder();

	void findConsensus();

};

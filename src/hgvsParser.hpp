#include <string>

class HgvsParser {
private:
	std::string & hgvsString;
	size_t posInUtr3;
public:
	HgvsParser(const std::string & hgvs);
	~HgvsParser();

private:
	void parse();
};

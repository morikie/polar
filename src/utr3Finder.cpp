#include <algorithm>
#include <string>
#include <vector>
#include <boost/foreach.hpp>
#include "utr3Finder.hpp"


/**
* Vector of known Poly(A) cleavage motifs to be not pathogenic. Sorted by frequency (highest first).
*/
const std::vector<std::string> Utr3Finder::hexamers = { 
	"aataaa",
	"attaaa",
	"tataaa",
	"agtaaa",
	"aagaaa",
	"aatata",
	"aataca",
	"cataaa",
	"gataaa",
	"aatgaa",
	"actaaa",
	"aataga"
//	"tttaaa",
//	"aaaaag",
//	"aaaaca",
//	"ggggct"
};


/**
* Vector with reversed hexamers. Same as hexamers, but reversed (not complemented). Used for searching with reverse_iterator.
*/
const std::vector<std::string> Utr3Finder::rHexamers( [] ()
{
	std::vector<std::string> temp(Utr3Finder::hexamers.begin(), Utr3Finder::hexamers.end());
	BOOST_FOREACH(std::string & str, temp) std::reverse(str.begin(), str.end());
	return temp;
}());


Utr3Finder::Utr3Finder(const SeqStruct & sSt):
	seqStruct(sSt)
{}


Utr3Finder::~Utr3Finder() {}


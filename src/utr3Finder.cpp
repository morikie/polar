#include <string>
#include <vector>
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
	"tttaaa",
	"actaaa",
	"aataga",
	"aaaaag",
	"aaaaca",
	"ggggct"
};


/**
* Vector with reversed hexamers. Same as hexamers, but reversed (not complemented). Used for searching with reverse_iterator.
*/
const std::vector<std::string> Utr3Finder::rHexamers = { 
	"aaataa",
	"aaatta",
	"aaatat",
	"aaatga",
	"aaagaa",
	"atataa",
	"acataa",
	"aaatac",
	"aaatag",
	"aagtaa",
	"aaattt",
	"aaatca",
	"agataa",
	"gaaaaa",
	"acaaaa",
	"tcgggg"
};


Utr3Finder::Utr3Finder(const SeqStruct & sSt):
	seqStruct(sSt)
{}


Utr3Finder::~Utr3Finder() {}


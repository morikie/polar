#include <climits>
#include <fstream>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/filesystem.hpp>
#include <seqan/seq_io.h>
#include "refGeneParser.hpp"
#include "polarUtility.hpp"


namespace polar {
namespace utility {

namespace qi = boost::spirit::qi;
namespace fs = boost::filesystem;


/**
 * Maps PAS motives to integers.
 */
size_t motifToIndex(std::string & motif) {
	if (motif == "aataaa") return 0;
	else if (motif == "attaaa") return 1;
	else if (motif == "tataaa") return 2;
	else if (motif == "agtaaa") return 3;
	else if (motif == "aagaaa") return 4;
	else if (motif == "aatata") return 5;
	else if (motif == "aataca") return 6;
	else if (motif == "cataaa") return 7;
	else if (motif == "gataaa") return 8;
	else if (motif == "aatgaa") return 9;
	else if (motif == "actaaa") return 10;
	else if (motif == "aataga") return 11;
	else return 12;
}


/**
 * Returns the complementary DNA base.
 */
char complement(const char c) {
	switch(c) {
	case 'a':
		return 't';
		break;
	case 'c':
		return 'g';
		break;
	case 'g':
		return 'c';
		break;
	case 't':
		return 'a';
		break;
	default:
		return 'n';
		break;
	}
}


/**
 * Maps the words "chr1", "chr2" etc. to integers for use as identifiers in the FastaIndex read functions.
 */
size_t getFastaIndex(const std::string & chr) {
	size_t idx = UINT_MAX;
	if (chr == "chrX" || chr == "X" || chr == "x") {
		idx = 22;
	} else if (chr == "chrY" || chr == "Y" || chr == "y") {
		idx = 23;
	} else {
		if (! qi::parse(chr.begin(), chr.end(), qi::omit[*qi::alpha] >> qi::uint_, idx)) {
			return UINT_MAX; 
		}
		//adjusting idx to 0-starting map
		idx--;
	}
	return idx; 
}


/**
 * Crude function to obtain the patch version of RefSeq transcripts from file ucsc_txRefSeq.txt (not offical file).
 */
std::unordered_map<std::string, size_t> getTxRefSeqAccessions(const fs::path & f) {
	typedef std::string accession;
	typedef size_t patch;
	std::unordered_map<accession, patch> txRefSeq;
	std::fstream in(f.string());
	std::string line;
	qi::rule<std::string::iterator, std::pair<accession, patch>() > r = *~qi::char_('\t') >> '\t' >> qi::uint_;
	while (std::getline(in, line)) {
		std::pair<accession, patch> temp;
		qi::parse(line.begin(), line.end(), r, temp);
		if (txRefSeq.find(temp.first) == txRefSeq.end()) txRefSeq[temp.first] = temp.second;
		else std::cerr << "Different patch numbers used for the same transcript" << std::endl;
	}
	return txRefSeq;
}


/**
 * Return the length of a 3' UTR.
 */
size_t getUtrLength(const RefGeneProperties & txProp) {
	size_t utr3Length = 0;
	if (txProp.strand == "+") {

		for (int i = txProp.exonStarts.size() - 1; i >= 0; i--) {
			if (txProp.cdsEnd >= txProp.exonStarts[i]) {
				utr3Length = txProp.exonEnds[i] - txProp.cdsEnd;
				for (size_t j = static_cast<size_t>(i) + 1; j < txProp.exonStarts.size(); j++) { 
					utr3Length += txProp.exonEnds[j] - txProp.exonStarts[j];
				}
				break;
			}
		}
	} else {
		for (size_t i = 0; i < txProp.exonEnds.size(); i++) {
			if (txProp.cdsStart <= txProp.exonEnds[i]) {
				utr3Length = txProp.cdsStart - txProp.exonStarts[i];
				for (int j = static_cast<int>(i) - 1; j >= 0; j--) { 
					utr3Length += txProp.exonEnds[j] - txProp.exonStarts[j];
				}
				break;
			}
		}
	}
	return utr3Length;
}


/**
 * Return the 3' UTR sequence.
 */
std::string getUtrSequence(const RefGeneProperties & txProp, seqan::FaiIndex & fai) {
	std::string returnString = std::string();
	size_t chr = polar::utility::getFastaIndex(txProp.chr);
	const std::string & strand = txProp.strand;
	if (strand == "-") {
		for (size_t i = 0; i < txProp.exonEnds.size(); i++) {
			if (txProp.cdsStart <= txProp.exonEnds[i]) {
				seqan::CharString finalString;
				seqan::CharString temp;
				for (size_t j = 0; j <= i; j++) {
					if (j == i) {
						seqan::readRegion(temp, fai, chr, txProp.exonStarts[j], txProp.cdsStart);
						seqan::append(finalString, temp);
						seqan::clear(temp);
					} else {
						seqan::readRegion(temp, fai, chr, txProp.exonStarts[j], txProp.exonEnds[j]);
						seqan::append(finalString, temp);
						seqan::clear(temp);
					}
				}
				seqan::reverseComplement(finalString);
				seqan::toLower(finalString);
				returnString = std::string(seqan::toCString(finalString));
				break;
			}
		}
	} else {
		for (int i = txProp.exonStarts.size() - 1; i >= 0; i--) {
			if (txProp.cdsEnd >= txProp.exonStarts[i]) {
				seqan::CharString finalString;
				seqan::CharString temp;
				seqan::readRegion(finalString, fai, chr, txProp.cdsEnd, txProp.exonEnds[i]);
				for (size_t j = i + 1; j < txProp.exonStarts.size(); j++) {
					seqan::readRegion(temp, fai, chr, txProp.exonStarts[j], txProp.exonEnds[j]);
					seqan::append(finalString, temp);
					seqan::clear(temp);
				}
				returnString = std::string(seqan::toCString(finalString));
				break;
			}
			
		}
	}
	return returnString;
}

/**
 * Mapping the genomic position to the transcript position. Integrity is checked.
 */
size_t mapGenomePosToTxPos(const RefGeneProperties & txProp, const size_t genomePos) {
	if (genomePos < txProp.exonStarts.front() || genomePos > txProp.exonEnds.back()) return UINT_MAX;
	const std::string & strand = txProp.strand;
	size_t txLength = polar::utility::getTxLength(txProp);
	size_t txPos = 0;
	for (size_t i = 0; i < txProp.exonEnds.size(); i++) {
		//inside another exon
		if (txProp.exonEnds[i] <= genomePos) {
			txPos += txProp.exonEnds[i] - txProp.exonStarts[i];
		//inside this exon
		} else if (txProp.exonEnds[i] > genomePos && txProp.exonStarts[i] <= genomePos) {
			txPos += genomePos - txProp.exonStarts[i];
			break;
		//genomePos is part of an intron; returning error code
		} else {
			return UINT_MAX;
		}
	}
	if (strand == "-") {
		return  txLength - txPos;
	} else {
		return txPos;
	}	
}


/** 
 * Calculate transcript length by summing over the exon lengths.
 */
size_t getTxLength(const RefGeneProperties & txProp) {
	size_t txLength = 0;
	for (size_t i = 0; i < txProp.exonStarts.size(); i++) txLength += txProp.exonEnds[i] - txProp.exonStarts[i];
	return txLength;
}


} /* namespace utility */
} /* namespace polar */


#define BOOST_TEST_MODULE polar_test
#define BOOST_TEST_MAIN polar 

#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <boost/optional/optional_io.hpp>
#include <boost/test/included/unit_test.hpp>
#include "../src/fastaReader.hpp"
#include "../src/jannovarVcfParser.hpp"
#include "../src/knownGeneParser.hpp"
#include "../src/polar.hpp"
#include "../src/readTranscripts.hpp"
#include "../src/seqStruct.hpp"
#include "../src/utr3Finder.hpp"
#include "../src/utr3FinderFuzzy.hpp"
#include "../src/utr3FinderNaive.hpp"


namespace fs = boost::filesystem;
namespace spirit = boost::spirit;

/**
 * Floating point comparison. Returns true if two floats are equal up to the desired ULP (units in the last place). 
 * Copied from http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon;
 */
template<typename T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
	almost_equal(T x, T y, int ulp)
{
	return std::abs(x-y) < std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
		|| std::abs(x-y) < std::numeric_limits<T>::min();
}

/*
BOOST_AUTO_TEST_CASE( fastaReader ) {
	fs::path f = "UTR3_sequences_test.txt";
	FastaReader newReader(f);
	
	BOOST_CHECK_EQUAL(newReader.getFile().string(), "UTR3_sequences_test.txt");
	
	FastaReader::fastaVector::const_iterator itBegin = newReader.getBeginIterator();
	FastaReader::fastaVector::const_iterator itEnd = newReader.getEndIterator();
	BOOST_CHECK_EQUAL(itBegin->first, 
"ENSG00000000003|ENST00000373020|;100630759;100628670|;10063\
0797;100629986|tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11\
858]|100628670|100636806");
	BOOST_CHECK_EQUAL(itBegin->second, 
"CCCAATGTATCTGTGGGCCTATTCCTCTCTACCTTTAAGGACATTTAGGGTCCCCCCTGT\
GAATTAGAAAGTTGCTTGGCTGGAGAACTGACAACACTACTTACTGATAGACCAAAAAAC\
TACACCAGTAGGTTGATTCAATCAAGATGTATGTAGACCTAAAACTACACCAATAGGCTG\
ATTCAATCAAGATCCGTGCTCGCAGTGGGCTGATTCAATCAAGATGTATGTTTGCTATGT\
TCTAAGTCCACCTTCTATCCCATTCATGTTAGATCGTTGAAACCCTGTATCCCTCTGAAA\
CACTGGAAGAGCTAGTAAATTGTAAATGAAGTAATACTGTGTTCCTCTTGACTGTTATTT\
TTCTTAGTAGGGGGCCTTTGGAAGGCACTGTGAATTTGCTATTTTGATGTAGTGTTACAA\
GATGGAAAATTGATTCCTCTGACTTTGCTATTGATGTAGTGTGATAGAAAATTCACCCCT\
CTGAACTGGCTCCTTCCCAGTCAAGGTTATCTGGTTTGATTGTATAATTTGCACCAAGAA\
GTTAAAATGTTTTATGACTCTCTGTTCTGCTGACAGGCAGAGAGTCACATTGTGTAATTT\
AATTTCAGTCAGTCAATAGATGGCATCCCTCATCAGGGTTGCCAGATGGTGATAACAGTG\
TAAGGCCTTGGGTCTAAGGCATCCACGACTGGAAGGGACTACTGATGTTCTGTGATACAT\
CAGGTTTCAGCACACAACTTACATTTCTTTGCCTCCAAATTGAGGCATTTATTATGATGT\
TCATACTTTCCCTCTTGTTTGAAAGTTTCTAATTATTAAATGGTGTCGGAATTGTTGTAT\
TTTCCTTAGGAATTCAGTGGAACTTATCTTCATTAAATTTAGCTGGTACCAGGTTGATAT\
GACTTGTCAATATTATGGTCAACTTTAAGTCTTAGTTTTCGTTTGTGCCTTTGATTAATA\
AGTATAACTCTTATACAATAAATACTGCTTTCCTCTAAAAAGATCGTGTTTAAATTAACT\
TGTAGAAAATCTGCTGGAATGGTTGTTGTTTTCCACTGAGAAAGCTAAGCCCTACATTTC\
TATTCAGAGTACTGTTTTTAGATGTGAAATATAAGCCTGCGGCCTTAACTCTGTATTAAA\
AAAAATGTTTTTGTTTAAAAAAAACTGTTCCCATAGGTGCAGCAAACCACCATGGCACAT\
GTATACCTATGTAACAAACCTGCACATTCTGCACATGTATCCCAGAACTTAATGTAAACA\
AAAAAATCTTAAAGTGCAAATATTAAAAAAAACTGTTCTCTGTGAAAAAAATTATATTCC\
ATGTTATAAAGTAGCATATGACTAGTGTTCTCCTAG");
	itBegin++;
	BOOST_CHECK_EQUAL(itBegin->first, 
"ENSG00000000003|ENST00000494424|||tetraspanin 6 [Source:HGN\
C Symbol;Acc:HGNC:11858]|100633442|100639991");
	BOOST_CHECK_EQUAL(itBegin->second, "Sequenceunavailable");
	itBegin++;
 	BOOST_CHECK_EQUAL(itBegin->first, 
"ENSG00000000003|ENST00000496771|||tetraspanin 6 [Source:HGN\
C Symbol;Acc:HGNC:11858]|100632541|100636689");
 	BOOST_CHECK_EQUAL(itBegin->second, "Sequenceunavailable");
	itBegin++;
	BOOST_CHECK_EQUAL(itBegin->first, 
"ENSG00000000003|ENST00000612152|;100630759;100627109|;10063\
0797;100629986|tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11\
858]|100627109|100637104");
	BOOST_CHECK_EQUAL(itBegin->second, 
"CCCAATGTATCTGTGGGCCTATTCCTCTCTACCTTTAAGGACATTTAGGGTCCCCCCTGT\
GAATTAGAAAGTTGCTTGGCTGGAGAACTGACAACACTACTTACTGATAGACCAAAAAAC\
TACACCAGTAGGTTGATTCAATCAAGATGTATGTAGACCTAAAACTACACCAATAGGCTG\
ATTCAATCAAGATCCGTGCTCGCAGTGGGCTGATTCAATCAAGATGTATGTTTGCTATGT\
TCTAAGTCCACCTTCTATCCCATTCATGTTAGATCGTTGAAACCCTGTATCCCTCTGAAA\
CACTGGAAGAGCTAGTAAATTGTAAATGAAGTAATACTGTGTTCCTCTTGACTGTTATTT\
TTCTTAGTAGGGGGCCTTTGGAAGGCACTGTGAATTTGCTATTTTGATGTAGTGTTACAA\
GATGGAAAATTGATTCCTCTGACTTTGCTATTGATGTAGTGTGATAGAAAATTCACCCCT\
CTGAACTGGCTCCTTCCCAGTCAAGGTTATCTGGTTTGATTGTATAATTTGCACCAAGAA\
GTTAAAATGTTTTATGACTCTCTGTTCTGCTGACAGGCAGAGAGTCACATTGTGTAATTT\
AATTTCAGTCAGTCAATAGATGGCATCCCTCATCAGGGTTGCCAGATGGTGATAACAGTG\
TAAGGCCTTGGGTCTAAGGCATCCACGACTGGAAGGGACTACTGATGTTCTGTGATACAT\
CAGGTTTCAGCACACAACTTACATTTCTTTGCCTCCAAATTGAGGCATTTATTATGATGT\
TCATACTTTCCCTCTTGTTTGAAAGTTTCTAATTATTAAATGGTGTCGGAATTGTTGTAT\
TTTCCTTAGGAATTCAGTGGAACTTATCTTCATTAAATTTAGCTGGTACCAGGTTGATAT\
GACTTGTCAATATTATGGTCAACTTTAAGTCTTAGTTTTCGTTTGTGCCTTTGATTAATA\
AGTATAACTCTTATACAATAAATACTGCTTTCCTCTAAAAAGATCGTGTTTAAATTAACT\
TGTAGAAAATCTGCTGGAATGGTTGTTGTTTTCCACTGAGAAAGCTAAGCCCTACATTTC\
TATTCAGAGTACTGTTTTTAGATGTGAAATATAAGCCTGCGGCCTTAACTCTGTATTAAA\
AAAAATGTTTTTGTTTAAAAAAAACTGTTCCCATAGGTGCAGCAAACCACCATGGCACAT\
GTATACCTATGTAACAAACCTGCACATTCTGCACATGTATCCCAGAACTTAATGTAAACA\
AAAAAATCTTAAAGTGCAAATATTAAAAAAAACTGTTCTCTGTGAAAAAAATTATATTCC\
ATGTTATAAAGTAGCATATGACTAGTGTTCTCCTAGAGATCAGACTTTTTTGATTGTATA\
GTTTGCATTAAAAAGTTGTACAGGGAGGGATGTAACCTGTATCTTCAGGATAATAGGGAA\
ATTAATAAGGAAAATAATAATTACTAAAATTTGAGTTGAAGTCAGGGAATTATTTCTTTG\
GTTTTGAGTCTCTTTATACATCCATTAGTAGAACCTGTCTAGTCTGATTGCCACAGTCCT\
TGAATTGATGGTAAGGGGGAGTCAGTTACAGTTAGAAAAAAACCTGGGCAAAAACTACTA\
GTTAAATGATCATGAATTTTGAGTATGTGTTTTAAAATGCTTGGAGTTATTGCAGAAAAA\
TGATATTCTCCATTAGAACATTTGATGTGTGACTTTTGACATAGGATAAGGTGATAAGAA\
AATTAAGTTGAAAAATAGTGACTCAAGTACAGCAAATCAGATTTTTGGCTACTGAGATTA\
AAGAAGCAAATATTTAAAGGATTATCAGCAATTTTCAAATACATATTTTTCACTGCAGTG\
GTTTTGAAGAATGTTTTAAGGTTATATTTAGAGTTTTTAAAAAAGTTTTCCCTCAGCCAG\
CCATGGTGGCTCACATCTGTAATCCCAGCACTTTGAGAGGCTGAGGCGGGCGGATCACCT\
GAGGTCAGGAGTTCTAGACCAGCCTGGCCAACATGGTAAAACACCGTCTCTAATAAAAAT\
GCAAAAATTAGCTGGGTGTGGTGGCAAGCACCTGTAATCCCAGCTACTTGGGAGACTGAG\
GCAGGAGAATCACTTGAATCCAGGAGTCAGAGGTTGCAGTGAGCCGAGATCACACCACTG\
CACTCCAGCCTGGGCGATAAGAGTGAACCTGCATCTCAGAAAAAAAAAAAAAGTTTTCCC\
TCATTATTAAAAAGGAAAAATTATAGAAAATTAAAATATAGAATGTAGAAAGAGGAAAAA\
CACCACACATCACTCGTCAGGTGGATTGGTTTGTTTATTTTTAAGATGGGTTTATTTTGG\
ATAGCCAGTTAGAAAACGCCCTTATAGCTGATGTGTGCTTGTGCATACTGAATTGTATAC\
TTAATTCTGAAAGTAATGAGGAAACATTTCACCTAATGATAACTGGTACAAAGGAAAGTT\
CACGTGTCTTTAATAGTTTGCTAATTAGACTAGCTGGTTAGAAAGAGGGATCACTAGATT\
GGAGAGAGAGATCCAGTTTCCACTGGAGGCAGTAGGTCTGCTCATTTCCATGGTGACAAT\
CCTAATGGCAGAGGGAGTGAGTTTCCCTTAGTGGCACAGGCTGGGGGTTTTGAGCTTTTA\
AGTTTTGTTTTATAGAAATCTGGAGAATTTTACATTTTACTTCACTATATATCTATGTTT\
AAAATAATTCATTTTTGATGTCCCAAAATTGATGTGGGTTGCCTTATTTCCAGGCACCAG\
GTTGCTCACTGCAAGGTGACATTTTGAGAAGCAGAAATCTCTAGGTATATCTTAAGTGTG\
AAAAGCTTGTTAGATTTCCATGGCCTATTCCAGGTGGGTTGTTGGGTTTATGTTGTAAGC\
ATCAAACTTTTGTGGAAATAAATACCTAATACAAGGT");
	itBegin++;
	BOOST_CHECK_EQUAL(itBegin->first, "ENSG00000000003|E\
NST00000614008|||tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:\
11858]|100632063|100637104");
	BOOST_CHECK_EQUAL(itBegin->second, "Sequenceunavailable");
	itBegin++;
	BOOST_CHECK_EQUAL(itBegin->first, "ENSG00000000005|E\
NST00000373031|;100599718|;100599885|tenomodulin [Source:HGN\
C Symbol;Acc:HGNC:17757]|100584802|100599885");
	BOOST_CHECK_EQUAL(itBegin->second, 
"TAGGAGGTTTGAGCTCAAATGCTTAAACTGCTGGCAACATATAATAAATGCATGCTATTC\
AATGAATTTCTGCCTATGAGGCATCTGGCCCCTGGTAGCCAGCTCTCCAGAATTACTTGT\
AGGTAATTCCTCTCTTCATGTTCTAATAAACTTCTACATTATCACCAA");
	itBegin++;
	BOOST_CHECK_EQUAL(itBegin->first, 
"ENSG00000000005|ENST00000485971|||tenomodulin [Source:HGNC \
Symbol;Acc:HGNC:17757]|100593624|100597531");
	BOOST_CHECK_EQUAL(itBegin->second, "Sequenceunavailable");
	itBegin++;
	BOOST_CHECK_EQUAL(itBegin->first, 
"ENSG00000000419|ENST00000371582|;50934867|;50935131|dolichy\
l-phosphate mannosyltransferase polypeptide 1, catalytic sub\
unit [Source:HGNC Symbol;Acc:HGNC:3005]|50934867|50958555");
	BOOST_CHECK_EQUAL(itBegin->second, 
"AAGAAAGATACTCATTTATAGTTACGTTCATTTCAGGTTAAACATGAAAGAAGCCTGGTT\
ACTGATTTGTATAAAATGTACTCTTAAAGTATAAAATATAAGGTAAGGTAAATTTCATGC\
ATCTTTTTATGAAGACCACCTATTTTATATTTCAAATTAAATAATTTTAAAGTTGCTGGC\
CTAATGAGCAATGTTCTCAATTTTCGTTTTCATTTTGCTGTATTGAGACCTATAAATAAA\
TGTATATTTTTTTTTGCATAAAGTA");
	itBegin++;
	BOOST_CHECK(itBegin == itEnd);

	
}
*/


/**
 * Tests for the knownGeneMrna Parser.
 */
BOOST_AUTO_TEST_CASE( knownGeneMrnaParser ) {
	
	fs::path file("ucsc_data/knownGeneTxMrna.txt");
	ReadTranscripts transcripts(file);
	
	BOOST_CHECK_EQUAL(transcripts.getValueByKey("uc001aaa.3"), 
"cttgccgtcagccttttctttgacctcttctttctgttcatgtgtatttgctgtctctta\
gcccagacttcccgtgtcctttccaccgggcctttgagaggtcacagggtcttgatgctg\
tggtcttcatctgcaggtgtctgacttccagcaactgctggcctgtgccagggtgcaagc\
tgagcactggagtggagttttcctgtggagaggagccatgcctagagtgggatgggccat\
tgttcatcttctggcccctgttgtctgcatgtaacttaataccacaaccaggcatagggg\
aaagattggaggaaagatgagtgagagcatcaacttctctcacaacctaggccagtgtgt\
ggtgatgccaggcatgcccttccccagcatcaggtctccagagctgcagaagacgacggc\
cgacttggatcacactcttgtgagtgtccccagtgttgcagaggcagggccatcaggcac\
caaagggattctgccagcatagtgctcctggaccagtgatacacccggcaccctgtcctg\
gacacgctgttggcctggatctgagccctggtggaggtcaaagccacctttggttctgcc\
attgctgctgtgtggaagttcactcctgccttttcctttccctagagcctccaccacccc\
gagatcacatttctcactgccttttgtctgcccagtttcaccagaagtaggcctcttcct\
gacaggcagctgcaccactgcctggcgctgtgcccttcctttgctctgcccgctggagac\
ggtgtttgtcatgggcctggtctgcagggatcctgctacaaaggtgaaacccaggagagt\
gtggagtccagagtgttgccaggacccaggcacaggcattagtgcccgttggagaaaaca\
ggggaatcccgaagaaatggtgggtcctggccatccgtgagatcttcccagggcagctcc\
cctctgtggaatccaatctgtcttccatcctgcgtggccgagggccaggcttctcactgg\
gcctctgcaggaggctgccatttgtcctgcccaccttcttagaagcgagacggagcagac\
ccatctgctactgccctttctataataactaaagttagctgccctggactattcaccccc\
tagtctcaatttaagaagatccccatggccacagggcccctgcctgggggcttgtcacct\
cccccaccttcttcctgagtcattcctgcagccttgctccctaacctgccccacagcctt\
gcctggatttctatctccctggcttggtgccagttcctccaagtcgatggcacctccctc\
cctctcaaccacttgagcaaactccaagacatcttctaccccaacaccagcaattgtgcc\
aagggccattaggctctcagcatgactatttttagagaccccgtgtctgtcactgaaacc\
ttttttgtgggagactattcctcccatctgcaacagctgcccctgctgactgcccttctc\
tcctccctctcatcccagagaaacaggtcagctgggagcttctgcccccactgcctaggg\
accaacaggggcaggaggcagtcactgaccccgagacgtttgcatcctgcacagctagag\
atcctttattaaaagcacactgttggtttctg");


	BOOST_CHECK_EQUAL(transcripts.getValueByKey("uc021ogw.1"), 
"gcgttggtggtttagtggtagaattctcgcctcccatgcgggagacccgggttcaattcc\
cggccactgca");

	BOOST_CHECK_EQUAL(transcripts.getValueByKey("uc011ncc.1"),
"cttgccgtcagccttttctttgacctcttctttctgttcatgtgtatttgctgtctctta\
gcccagacttcccgtgtcctttccaccaggcctttgagaggtcacagggtcttgatgctgt\
ggtcttgatctgcaggtgtctgacttccagcaactgctggcctgtgccagggtgcaagctg\
agcactggagtggagttttcctgtggagaggagccatgcctagagtgggatgggccattgt\
tcatcttctggcccctgttgtctgcatgtaacttaataccacaaccaggcataggggaaag\
attggaggaaagatgagtgagagcatcaacttctctgacaacctaggccagtgtgtggtga\
tgccaggcatgcccttccccagcatcaggtctccagagctgcagaagacgacggccgactt\
ggatcacaatcttgtgagtgtccccagtgttgcagaggcagggccatcaggcaccaaaggg\
attctgccagcatagtgctcctggattagtgatacacccggcaccctgtcctggacaagct\
gttggcctggatctgagccctcgtggaggtcaaagccacctttggttctgccattgctgct\
gtgtggaagttcactcctgccttttcctttccctagagcctccaccaccccgagatcacat\
ttctcactgccttttgtctgcccagtttcaccagaagtaggcctcttcctgacaggcagct\
gcaccactgcctggcgctgcgcccttcctttgctctgcccgctggagacggtgtttgtcat\
gggcctgatctgcagggatcctgctacaaaggtgaaacccagaagagtgtggagtccagag\
tgttgccaggacccaggcacaggcattagtgcccgttggagaaaacaggggaaccccgaag\
aaatggtgggtcctggccatccgtgagatcttcccagggcagctcccctctgtggaatcca\
atctgtcttccatcctgtgtggccgagggccaggcttctcactgggcctctgcaggaggct\
gccatttgtcctgcccaccttcttagaagcgagacggagcagacccatctgctactgccct\
ttctataataactaaagttagctgccctggactattcaccccctagtctcaatttaaaaag\
atccccatggccacagggcccctgcctgggggcttgtcacctcccccaccttcttcctgag\
tcacccctgcagccttgctccctaacctgccccacagccttgcctggatttctatctccct\
ggcttggtgccagttcctccaagttgatggcacctccctccctctcaaccacttgagcaaa\
ctccaagacatcttctaccccaacaccagcaattgtgccaagggccattgggctctcagca\
tgactatttttagagaccccgtgtctgtcactgaaaccttttttgtgggagactattcctc\
ccatctgcaacagctgcccctgctgactgcccttctctcccagagaaacaggtcagctggg\
agcttctgcccccactgcctagggaccaacaggggcaggaggcagtcactgaccccgagac\
gtttgcatcctgcacagctagaggtcctttattaaaagcacactgttggtttctgctc");
}


/**
 * Tests for the knownGene Parser.
 */
BOOST_AUTO_TEST_CASE( knownGeneParser ) {

	fs::path file = "ucsc_data/knownGene.txt";

	KnownGeneParser newParser(file);
	
	TxProperties txP1 { 
		"chr1",
		"+",
		11873,
		14409,
		11873,
		11873,
		std::vector<unsigned int> {11873, 12612, 13220},
		std::vector<unsigned int> {12227,12721,14409}
	};
	BOOST_CHECK(newParser.getValueByKey("uc001aaa.3") == txP1); 

	TxProperties txP2 {
		"chr12",
		"+",
		72666528,
		73059422,
		72666558,
		73056975,
		std::vector<unsigned int> {72666528,72680460,72771774,72863537,72866846,72893277,72936070,72955944,
			72956632,72962347,72969034,72969266,73012670,73014887,73015423,73046101,73046795,73050706,73056831},
		std::vector<unsigned int> {72667337,72680734,72771901,72863692,72866960,72893415,72936136,72956010,
			72956820,72962436,72969168,72969322,73012818,73014985,73015531,73046269,73046936,73050788,73059422}
	};
	BOOST_CHECK(newParser.getValueByKey("uc001sxa.3") == txP2);

	TxProperties txP3 {
		"chrY",
		"-",
		59358328,
		59360854,
		59358328,
		59358328,
		std::vector<unsigned int> {59358328,59360006,59360500},
		std::vector<unsigned int> {59359508,59360115,59360854}
	};
	BOOST_CHECK(newParser.getValueByKey("uc011ncc.1") == txP3);
}


/**
 * Tests for the Jannovar VCF Parser.
 */
BOOST_AUTO_TEST_CASE( jannovarVcfParser ) {
	fs::path file = "vcf/vcf-example.jv.vcf";
	JannovarVcfParser newParser(file);

	vcfTranscripts vcfTx1 { 
		"3_prime_utr_variant",
		"uc022cgy.1",
		"c.*197_*201del"
	};
	auto key1 = std::make_pair("X", 151283436u);
	BOOST_CHECK(newParser.getData().find(key1)->second[0] == vcfTx1);
	
	vcfTranscripts vcfTx1_1 {
		"3_prime_utr_variant",
		"uc004ffj.3",
		"c.*197_*201del"
	};
	auto key1_1 = std::make_pair("X", 151283436u);
	BOOST_CHECK(newParser.getData().find(key1_1)->second[1] == vcfTx1_1);

	vcfTranscripts vcfTx2 { 
		"non_coding_transcript_intron_variant",
		"uc021wmm.1",
		"n.107+55407T>C"
	};
	auto key2 = std::make_pair("22", 22842206u);
	BOOST_CHECK(newParser.getData().find(key2)->second[2] == vcfTx2);
	
	vcfTranscripts vcfTx3 { 
		"missense_variant",
		"uc031pjn.1",
		"c.508G>C"
	};
	auto key3 = std::make_pair("1", 879482u);
	BOOST_CHECK(newParser.getData().find(key3)->second[0] == vcfTx3);
	
	vcfTranscripts vcfTx4 { 
		"",
		"",
		""
	};
	auto key4 = std::make_pair("GL000229.1", 13276u);
	BOOST_CHECK(newParser.getData().find(key4)->second[0] == vcfTx4);
	
	vcfTranscripts vcfTx5 {
		"non_coding_transcript_exon_variant",
		"uc031pkn.1",
		"n.1727C>T"
	};
	auto key5 = std::make_pair("1", 879317u);
	BOOST_CHECK (newParser.getData().find(key5)->second[32] == vcfTx5);
}


/**
 * Tests for the HGVS string parser. 
 */
BOOST_AUTO_TEST_CASE ( hgvsParser ) {
	std::string hgvs = "c.*150A>G";
	HgvsParser newHgvs1 (hgvs);
	BOOST_CHECK_EQUAL(newHgvs1.getMutPosition(), 149);
	BOOST_CHECK_EQUAL(newHgvs1.isIntronic(), false);

	hgvs = "c.*114122_*114124del";
	HgvsParser newHgvs2 (hgvs);
	BOOST_CHECK_EQUAL(newHgvs2.getMutPosition(), 114121);
	BOOST_CHECK_EQUAL(newHgvs2.isIntronic(), false);
	
	hgvs = "c.*150+20A>G";
	HgvsParser newHgvs3 (hgvs);
	BOOST_CHECK_EQUAL(newHgvs3.getMutPosition(), 149);
	BOOST_CHECK_EQUAL(newHgvs3.isIntronic(), true);

	hgvs = "c.*150-20A>G";
	HgvsParser newHgvs4 (hgvs);
	BOOST_CHECK_EQUAL(newHgvs4.getMutPosition(), 149);
	BOOST_CHECK_EQUAL(newHgvs4.isIntronic(), true);	
}


/**
 * Tests for class Utr3FinderNaive.
 */
BOOST_AUTO_TEST_CASE ( utr3FinderNaive ) {	
	std::string seq ("acaataaaacccccccccccccatttttttttttggggtagagatagagccgagcagata\
gcccagagcacagtataccaagagagaataaaccaaaaaaaaaaaaaaaaaaaa");
	boost::optional<const HgvsParser> newHgvs1 = HgvsParser("c.*68A>G");
	size_t utr3Start = 20;
	size_t txLength = seq.size() - 15;
	SeqStruct  txTest1 {
		seq,
		utr3Start,
		txLength,
		newHgvs1,
		boost::none,
		boost::none,
		boost::none,
		boost::none
	};
	size_t utr3MotifPos = 86;
	Utr3FinderNaive utr3MutFi_test1 (txTest1);
	BOOST_CHECK_EQUAL(utr3MutFi_test1.getPolyaMotifPos()[0].pos, utr3MotifPos);
	BOOST_CHECK_EQUAL(utr3MutFi_test1.getMotifSequence(utr3MutFi_test1.getPolyaMotifPos()[0]), std::string ("aataaa"));	
	BOOST_CHECK_EQUAL(utr3MutFi_test1.isMutationInMotif(), true);		

	seq = "acaaataatataccaagagagaataaaccaaaaaaaaaaaaaaaaaa";
	boost::optional<const HgvsParser> newHgvs2 = HgvsParser("c.*2A>G");
	utr3Start = 20;
	txLength = seq.size();
	SeqStruct  txTest2 {
		seq,
		utr3Start,
		txLength,
		newHgvs2,
		boost::none,
		boost::none,
		boost::none,
		boost::none
	};
	utr3MotifPos = 21;
	Utr3FinderNaive utr3MutFi_test2 (txTest2);
	BOOST_CHECK_EQUAL(utr3MutFi_test2.getPolyaMotifPos()[0].pos, utr3MotifPos);
	BOOST_CHECK_EQUAL(utr3MutFi_test2.getMotifSequence(utr3MutFi_test2.getPolyaMotifPos()[0]), std::string ("aataaa"));
	BOOST_CHECK_EQUAL(utr3MutFi_test2.isMutationInMotif(), true);	

	seq = "acaataaaacccccccccccccatttttttttttggggtagagatagagccgagcagatagcccagagcacagtatataaaccaagagagaaaaacc";
	boost::optional<const HgvsParser> newHgvs3 = HgvsParser("c.*21A>G");
	utr3Start = 60;
	txLength = seq.size();
	SeqStruct  txTest3 {
		seq,
		utr3Start,
		txLength,
		newHgvs3,
		boost::none,
		boost::none,
		boost::none,
		boost::none
	};
	utr3MotifPos = 75;
	Utr3FinderNaive utr3MutFi_test3 (txTest3);
	BOOST_CHECK_EQUAL(utr3MutFi_test3.getPolyaMotifPos()[0].pos, utr3MotifPos);
	BOOST_CHECK_EQUAL(utr3MutFi_test3.getMotifSequence(utr3MutFi_test3.getPolyaMotifPos()[0]), std::string ("tataaa"));
	BOOST_CHECK_EQUAL(utr3MutFi_test3.isMutationInMotif(), true);	

	seq = "accccaaatatccccccccacagtatataaaccaagagagaaaaacc";
	boost::optional<const HgvsParser> newHgvs4 = HgvsParser("c.*6A>G");
	utr3Start = 20;
	txLength = seq.size();
	SeqStruct  txTest4 {
		seq,
		utr3Start,
		txLength,
		newHgvs4,
		boost::none,
		boost::none,
		boost::none,
		boost::none
	};
	utr3MotifPos = 25;
	Utr3FinderNaive utr3MutFi_test4 (txTest4);
	BOOST_CHECK_EQUAL(utr3MutFi_test4.getPolyaMotifPos()[0].pos, utr3MotifPos);
	BOOST_CHECK_EQUAL(utr3MutFi_test4.getMotifSequence(utr3MutFi_test4.getPolyaMotifPos()[0]), std::string ("tataaa"));
	BOOST_CHECK_EQUAL(utr3MutFi_test4.isMutationInMotif(), true);	

	seq = "accccaatccccccccacagtataaccaagagagaaaaacc";
	boost::optional<const HgvsParser> newHgvs5 = HgvsParser("c.*2A>G");
	utr3Start = 10;
	txLength = seq.size();
	SeqStruct  txTest5 {
		seq,
		utr3Start,
		txLength,
		newHgvs5,
		boost::none,
		boost::none,
		boost::none,
		boost::none
	};
	utr3MotifPos = Utr3FinderNaive::noHitPos;
	Utr3FinderNaive utr3MutFi_test5 (txTest5);
	BOOST_CHECK_EQUAL (utr3MutFi_test5.getPolyaMotifPos()[0].pos, utr3MotifPos);
	BOOST_CHECK (utr3MutFi_test5.getMotifSequence(utr3MutFi_test5.getPolyaMotifPos()[0]) == "");
	BOOST_CHECK_EQUAL (utr3MutFi_test5.isMutationInMotif(), false);	
}


/* Tests for Utr3FinderFuzzy */

BOOST_AUTO_TEST_CASE ( utr3FinderFuzzy ) {	
	std::string seq = "agggagagattaatttttcccaccccaataaaccccccccacagtataaccaagagattttatttgaaaaacc";
	SeqStruct  txTest1 {
		seq,
		boost::none,
		boost::none,
		boost::none,
		boost::none,
		boost::none,
		boost::none,
		boost::none
	};
	size_t utr3MotifPos = 26;
	Utr3FinderFuzzy utr3FinderFuz_test1 (txTest1);
	auto pairLeft1 = Utr3FinderFuzzy::dseLocMap.find(std::string("aataaa"))->second.getLeftStraight();
	auto pairRight1 = Utr3FinderFuzzy::dseLocMap.find(std::string("aataaa"))->second.getRightStraight();
	double slopeLeft1 = 1.0 / 15;
	double interceptLeft1 = - 2.0 / 3;
	double slopeRight1 = - 1.0 / 20;
	double interceptRight1 = 11.0 / 4;
	std::cerr << pairLeft1.first << ", " << slopeLeft1 << std::endl;
	BOOST_CHECK(almost_equal(pairLeft1.first, slopeLeft1, 2));
	BOOST_CHECK(almost_equal(pairLeft1.second, interceptLeft1, 2));
	BOOST_CHECK(almost_equal(pairRight1.first, slopeRight1, 2));
	BOOST_CHECK(almost_equal(pairRight1.second, interceptRight1, 2));

	auto pairLeft2 = Utr3FinderFuzzy::dseLocMap.find(std::string("attaaa"))->second.getLeftStraight();
	auto pairRight2 = Utr3FinderFuzzy::dseLocMap.find(std::string("attaaa"))->second.getRightStraight();
	double slopeLeft2 = 1.0 / 15;
	double interceptLeft2 = - 2.0 / 3;
	double slopeRight2 = - 1.0 / 27;
	double interceptRight2 = 20.0 / 9;
	BOOST_CHECK(almost_equal(pairLeft2.first, slopeLeft2, 2));
	BOOST_CHECK(almost_equal(pairLeft2.second, interceptLeft2, 2));
	BOOST_CHECK(almost_equal(pairRight2.first, slopeRight2, 2));
	BOOST_CHECK(almost_equal(pairRight2.second, interceptRight2, 2));
	BOOST_CHECK_EQUAL(utr3MotifPos, utr3FinderFuz_test1.getPolyaMotifPos()[0].pos);
//	BOOST_CHECK_EQUAL (utr3FinderFuz_test1.getPolyaMotifPos()[0], utr3MotifPos);
//	BOOST_CHECK (utr3FinderFuz_test1.getMotifSequence(utr3FinderFuz_test1.getPolyaMotifPos()[0]) == "aataaa");
//	BOOST_CHECK_EQUAL (utr3FinderFuz_test1.isMutationInMotif(), true);	
}

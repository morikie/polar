#define BOOST_TEST_MAIN polar 

#include <iostream>
#include <unordered_map>
#include <boost/test/unit_test.hpp>
#include "../src/fastaReader.hpp"
#include "../src/jannovarVcfParser.hpp"
#include "../src/knownGeneParser.hpp"
#include "../src/polar.hpp"
#include "../src/readTranscripts.hpp"
#include "../src/transcriptMutation.hpp"
#include "../src/utr3MutationFinder.hpp"

namespace fs = boost::filesystem;
namespace spirit = boost::spirit;

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

BOOST_AUTO_TEST_CASE( knownGeneMrnaParser ) {
	
	fs::path file("knownGeneMrna.txt");
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


BOOST_AUTO_TEST_CASE( knownGeneParser ) {

	fs::path file = "knownGene.txt";

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


BOOST_AUTO_TEST_CASE( jannovarVcfParser ) {
	fs::path file = "vcf-example.jv.vcf";

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


BOOST_AUTO_TEST_CASE( utr3MutationFinder ) {	
	
	std::string chrom ("chr1");
	size_t gPos = 12345;
	std::string strand ("+");
	std::string seqId ("uc002wel.2");
	std::string seq ("acaataaaacccccccccccccatttttttttttggggtagagatagagccgagcagatagcccagagcacagtataccaagagagaataaacc");
	HgvsParser newHgvs1 ("c.*68A>G");
	size_t utr3Start = 20;
	TranscriptMutation  txTest1 {
		chrom,
		gPos,
		strand,
		seqId,
		seq,
		newHgvs1,
		utr3Start
	};
	size_t utr3MotifPos = 86;
	Utr3MutationFinder utr3MutFi_test1 (txTest1);
	BOOST_CHECK_EQUAL(utr3MutFi_test1.getPolyaMotifPos(), utr3MotifPos);
	BOOST_CHECK_EQUAL(utr3MutFi_test1.getMotifSequence(), std::string ("aataaa"));	
	BOOST_CHECK_EQUAL(utr3MutFi_test1.isMutationInMotif(), true);	
	

	chrom = "chr1";
	gPos = 12345;
	strand = "-";
	seqId = "uc002wel.2";
	seq = "acaaataatataccaagagagaataaacc";
	HgvsParser newHgvs2 ("c.*2A>G");
	utr3Start = 20;
	TranscriptMutation  txTest2 {
		chrom,
		gPos,
		strand,
		seqId,
		seq,
		newHgvs2,
		utr3Start
	};
	utr3MotifPos = 21;
	Utr3MutationFinder utr3MutFi_test2 (txTest2);
	BOOST_CHECK_EQUAL(utr3MutFi_test2.getPolyaMotifPos(), utr3MotifPos);
	BOOST_CHECK_EQUAL(utr3MutFi_test2.getMotifSequence(), std::string ("aataaa"));
	BOOST_CHECK_EQUAL(utr3MutFi_test2.isMutationInMotif(), true);	

	chrom = "chr1";
	gPos = 12345;
	strand = "+";
	seqId = "uc002wel.2";
	seq = "acaataaaacccccccccccccatttttttttttggggtagagatagagccgagcagatagcccagagcacagtatataaaccaagagagaaaaacc";
	HgvsParser newHgvs3 ("c.*21A>G");
	utr3Start = 60;
	TranscriptMutation  txTest3 {
		chrom,
		gPos,
		strand,
		seqId,
		seq,
		newHgvs3,
		utr3Start
	};
	utr3MotifPos = 75;
	Utr3MutationFinder utr3MutFi_test3 (txTest3);
	BOOST_CHECK_EQUAL(utr3MutFi_test3.getPolyaMotifPos(), utr3MotifPos);
	BOOST_CHECK_EQUAL(utr3MutFi_test3.getMotifSequence(), std::string ("tataaa"));
	BOOST_CHECK_EQUAL(utr3MutFi_test3.isMutationInMotif(), true);	

	chrom = "chr1";
	gPos = 12345;
	strand = "-";
	seqId = "uc002wel.2";
	seq = "accccaaatatccccccccacagtatataaaccaagagagaaaaacc";
	HgvsParser newHgvs4 ("c.*6A>G");
	utr3Start = 20;
	TranscriptMutation  txTest4 {
		chrom,
		gPos,
		strand,
		seqId,
		seq,
		newHgvs4,
		utr3Start
	};
	utr3MotifPos = 25;
	Utr3MutationFinder utr3MutFi_test4 (txTest4);
	BOOST_CHECK_EQUAL(utr3MutFi_test4.getPolyaMotifPos(), utr3MotifPos);
	BOOST_CHECK_EQUAL(utr3MutFi_test4.getMotifSequence(), std::string ("tataaa"));
	BOOST_CHECK_EQUAL(utr3MutFi_test4.isMutationInMotif(), true);	

	chrom = "chr1";
	gPos = 12345;
	strand = "-";
	seqId = "uc002wel.2";
	seq = "accccaatccccccccacagtataaccaagagagaaaaacc";
	HgvsParser newHgvs5 ("c.*2A>G");
	utr3Start = 10;
	TranscriptMutation  txTest5 {
		chrom,
		gPos,
		strand,
		seqId,
		seq,
		newHgvs5,
		utr3Start
	};
	utr3MotifPos = Utr3MutationFinder::noHitPos;
	Utr3MutationFinder utr3MutFi_test5 (txTest5);
	BOOST_CHECK_EQUAL (utr3MutFi_test5.getPolyaMotifPos(), utr3MotifPos);
	BOOST_CHECK (utr3MutFi_test5.getMotifSequence() == "");
	BOOST_CHECK_EQUAL (utr3MutFi_test5.isMutationInMotif(), false);	
	

}

#define BOOST_TEST_MAIN polar 

#include <iostream>
#include <unordered_map>
#include <boost/spirit/home/support/multi_pass.hpp>
#include <boost/test/unit_test.hpp>
#include "../src/fastaReader.hpp"
#include "../src/jannovarVcfParser.hpp"
#include "../src/knownGeneParser.hpp"
#include "../src/polar.hpp"
#include "../src/readTranscripts.hpp"
#include "../src/utr3MutationFinder.hpp"

namespace fs = boost::filesystem;
namespace spirit = boost::spirit;

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

BOOST_AUTO_TEST_CASE( knownGeneMrnaParser ) {
	
	fs::path file("knownGeneMrna.txt");
	ReadTranscripts transcripts(file);
	
	BOOST_CHECK_EQUAL(transcripts.getData()["uc001aaa.3"], 
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


	BOOST_CHECK_EQUAL(transcripts.getData()["uc021ogw.1"], 
"gcgttggtggtttagtggtagaattctcgcctcccatgcgggagacccgggttcaattcc\
cggccactgca");

	BOOST_CHECK_EQUAL(transcripts.getData()["uc011ncc.1"],
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
	
	TxProperties txP1{ 
		"chr1",
		"+",
		11873,
		14409,
		11873,
		11873,
		std::vector<unsigned int> {11873, 12612, 13220},
		std::vector<unsigned int> {12227,12721,14409}
	};


	BOOST_CHECK(newParser.getData()["uc001aaa.3"] == txP1); 
	
	std::cerr << "newParser.getData().size(): " << newParser.getData().size() << std::endl;
}

BOOST_AUTO_TEST_CASE( jannovarVcfParser ) {
	
	fs::path file = "vcf-example.jv.vcf";

	JannovarVcfParser newParser(file);

	vcfTranscripts vcfTx1 { 
		"coding_transcript_intron_variant",
		"uc002wel.4",
		"c.365+2545A>G"
	};
	BOOST_CHECK(newParser.getData()[std::make_pair("20", 1110696u)][0] == vcfTx1);

	
	vcfTranscripts vcfTx2 { 
		"intergenic_variant",
		"uc002wcw.3",
		""
	};
	BOOST_CHECK(newParser.getData()[std::make_pair("20", 14370u)][0] == vcfTx2);
	
	vcfTranscripts vcfTx3 { 
		"3_prime_utr_variant",
		"uc002wep.4",
		"c.*2682A>G"
	};
	BOOST_CHECK(newParser.getData()[std::make_pair("20", 1148406u)][2] == vcfTx3);
	
	vcfTranscripts vcfTx4 { 
		"",
		"",
		""
	};
	BOOST_CHECK(newParser.getData()[std::make_pair("20", 1230237u)][0] == vcfTx4);

}


BOOST_AUTO_TEST_CASE( findConsensus ) {	

}

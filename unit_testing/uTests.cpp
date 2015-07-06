#define BOOST_TEST_MAIN polar 

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <seqan/find.h>
#include "../src/polar.hpp"
#include "../src/fastaReader.hpp"
#include "../src/mutationFinder.hpp"

namespace fs = boost::filesystem;


BOOST_AUTO_TEST_CASE( readIn ) {
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


BOOST_AUTO_TEST_CASE( findConsensus ) {	
	fs::path f = "UTR3_sequences_test.txt";
	FastaReader newReader(f);
	
	FastaReader::fastaVector::const_iterator itBegin = newReader.getBeginIterator();
	FastaReader::fastaVector::const_iterator itEnd = newReader.getEndIterator();
	std::cerr << "in findConsensus" << std::endl;	
	MutationFinder muFi(itBegin->second);

}

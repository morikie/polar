CC 		= g++
CFLAGS		= -c -Wall -Wextra -pedantic -O2 -std=c++14 -g
LFLAGS		= -Wall -fopenmp
BOOST		= -I $(HOME)/boost_1_62_0 
SEQAN		= -I $(HOME)/seqan/include
VIENNA		= -I $(HOME)/Vienna/include
LIBPATH		= -Llib/
CXX		= /usr/bin/gcc
LIBS		= -lboost_regex -lboost_filesystem -lboost_system -lboost_iostreams -lRNA
LIBSTEST	= -lboost_unit_test_framework -lboost_filesystem -lboost_system -lboost_iostreams
TPATH		= bin/

OBJS		= $(TPATH)buildIndexFile.o \
	$(TPATH)buildSeqStruct.o \
	$(TPATH)fastaReader.o \
	$(TPATH)hgvsParser.o \
	$(TPATH)jannovarVcfParser.o \
	$(TPATH)knownGeneParser.o \
	lib/polarUtility.o \
	$(TPATH)readTranscripts.o \
	$(TPATH)seqStruct.o \
	$(TPATH)utr3Finder.o \
	$(TPATH)utr3FinderFuzzy.o \
	$(TPATH)utr3FinderNaive.o \

.PHONY : all
all : $(TPATH)polar $(TPATH)uTests $(TPATH)acc_test $(TPATH)perf_fuzzy $(TPATH)perf_bppFuzzy

lib/polarUtility.o : src/polarUtility.hpp src/polarUtility.cpp 
	@echo "[Compile] polar utility functions"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) src/polarUtility.cpp -o lib/polarUtility.o

$(TPATH)polar : $(TPATH)polar.o $(OBJS) 
	@echo "[Link] polar"
	@$(CC) $(BOOST) $(SEQAN) $(VIENNA) $^ $(LIBPATH) $(LFLAGS) $(LIBS) -o $(TPATH)polar 

#unit testing
$(TPATH)uTests : $(TPATH)uTests.o $(OBJS) 
	@echo "[Link] uTests"
	@$(CC) $^ $(LIBPATH) $(LFLAGS) $(LIBSTEST) -o bin/uTests

$(TPATH)uTests.o : unit_testing/uTests.cpp
	@echo "[Compile] uTests"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) unit_testing/uTests.cpp -o $(TPATH)uTests.o

#polar main program
$(TPATH)polar.o : src/polar.hpp \
	src/polar.cpp \
	src/buildSeqStruct.hpp \
	src/jannovarVcfParser.hpp \
	src/knownGeneParser.hpp \
	src/polarUtility.hpp \
	src/readTranscripts.hpp \
	src/seqStruct.hpp \
	src/utr3Finder.hpp \
	src/utr3FinderFuzzy.hpp \
	src/utr3FinderNaive.hpp 
	@echo "[Compile] polar"
	@$(CC) $(BOOST) $(SEQAN) $(VIENNA) $(CFLAGS) src/polar.cpp -o $(TPATH)polar.o

$(TPATH)buildIndexFile.o : src/buildIndexFile.hpp src/buildIndexFile.cpp
	@echo "[Compile] BuildIndexFile"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/buildIndexFile.cpp -o $(TPATH)buildIndexFile.o

$(TPATH)fastaReader.o : src/fastaReader.hpp src/fastaReader.cpp
	@echo "[Compile] FastaReader"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/fastaReader.cpp -o $(TPATH)fastaReader.o

$(TPATH)hgvsParser.o : src/hgvsParser.hpp src/hgvsParser.cpp
	@echo "[Compile] HGVS Parser"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/hgvsParser.cpp -o $(TPATH)hgvsParser.o

$(TPATH)jannovarVcfParser.o : src/jannovarVcfParser.hpp src/jannovarVcfParser.cpp
	@echo "[Compile] JannovarVcfParser"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/jannovarVcfParser.cpp -o $(TPATH)jannovarVcfParser.o

$(TPATH)knownGeneParser.o : src/knownGeneParser.hpp src/knownGeneParser.cpp
	@echo "[Compile] KnownGeneParser"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/knownGeneParser.cpp -o $(TPATH)knownGeneParser.o
$(TPATH)buildSeqStruct.o : src/buildSeqStruct.hpp src/buildSeqStruct.cpp
	@echo "[Compile] buildSeqStruct"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/buildSeqStruct.cpp -o $(TPATH)buildSeqStruct.o

$(TPATH)readTranscripts.o : src/readTranscripts.hpp src/readTranscripts.cpp
	@echo "[Compile] ReadTranscripts"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/readTranscripts.cpp -o $(TPATH)readTranscripts.o

$(TPATH)seqStruct.o : src/seqStruct.hpp src/seqStruct.cpp src/hgvsParser.hpp
	@echo "[Compile] SeqStruct"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/seqStruct.cpp -o $(TPATH)seqStruct.o

$(TPATH)utr3Finder.o : src/utr3Finder.hpp src/utr3Finder.cpp src/hgvsParser.hpp src/seqStruct.hpp
	@echo "[Compile] UTR3Finder"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/utr3Finder.cpp -o $(TPATH)utr3Finder.o

$(TPATH)utr3FinderNaive.o : src/utr3FinderNaive.hpp src/utr3FinderNaive.cpp src/utr3Finder.hpp
	@echo "[Compile] UTR3FinderNaive"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/utr3FinderNaive.cpp -o $(TPATH)utr3FinderNaive.o

$(TPATH)utr3FinderFuzzy.o : src/utr3FinderFuzzy.hpp src/utr3FinderFuzzy.cpp src/utr3Finder.hpp
	@echo "[Compile] UTR3FinderFuzzy"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/utr3FinderFuzzy.cpp -o $(TPATH)utr3FinderFuzzy.o

$(TPATH)bppPredictFuzzy.o : src/bppPredictFuzzy.hpp src/bppPredictFuzzy.cpp src/utr3FinderFuzzy.hpp src/utr3Finder.hpp
	@echo "[Compile] BPPPredictFuzzy"
	@$(CC) $(BOOST) $(SEQAN) $(VIENNA) $(LIBPATH) $(CFLAGS) $(LIBS) src/bppPredictFuzzy.cpp -o $(TPATH)bppPredictFuzzy.o

#targets for the performance tests
$(TPATH)acc_test : $(TPATH)acc_test.o \
	$(TPATH)hgvsParser.o \
	lib/polarUtility.o \
	$(TPATH)readKnownPolyA.o \
	$(TPATH)utr3Finder.o \
	$(TPATH)utr3FinderFuzzy.o \
	$(TPATH)utr3FinderNaive.o \
	$(TPATH)seqStruct.o
	@echo "[Link] acc_test"
	@$(CC) $(BOOST) $(SEQAN) $^ $(LIBPATH) $(LFLAGS) $(LIBS) -o $(TPATH)acc_test

$(TPATH)createTPset.o : perf_testing/createTPset.hpp perf_testing/createTPset.cpp
	@echo "[Compile] TP data set"
	@$(CC) $(BOOST) $(SEQAN) $(VIENNA) $(LIBPATH) $(CFLAGS) $(LIBS) perf_testing/createTPset.cpp -o $(TPATH)createTPset.o

$(TPATH)createTNset.o : perf_testing/createTNset.hpp perf_testing/createTNset.cpp
	@echo "[Compile] TN data set"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) perf_testing/createTNset.cpp -o $(TPATH)createTNset.o

$(TPATH)acc_test.o : perf_testing/acc_test.hpp perf_testing/acc_test.cpp perf_testing/readKnownPolyA.hpp
	@echo "[Compile] acc_test"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) perf_testing/acc_test.cpp -o $(TPATH)acc_test.o

$(TPATH)readKnownPolyA.o : perf_testing/readKnownPolyA.hpp perf_testing/readKnownPolyA.cpp 
	@echo "[Compile] readKnownPolyA"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) perf_testing/readKnownPolyA.cpp -o $(TPATH)readKnownPolyA.o

$(TPATH)perf_fuzzy : $(TPATH)perf_fuzzy.o \
	$(TPATH)createTPset.o \
	$(TPATH)createTNset.o \
	$(TPATH)hgvsParser.o \
	lib/polarUtility.o \
	$(TPATH)readKnownPolyA.o \
	$(TPATH)refGeneParser.o \
	$(TPATH)utr3Finder.o \
	$(TPATH)utr3FinderFuzzy.o \
	$(TPATH)utr3FinderNaive.o \
	$(TPATH)seqStruct.o 
	@echo "[Link] perf_fuzzy"
	@$(CC) $(BOOST) $(SEQAN) $^ $(LIBPATH) $(LFLAGS) $(LIBS) -o $(TPATH)perf_fuzzy

$(TPATH)perf_fuzzy.o : perf_testing/perf_fuzzy.hpp perf_testing/perf_fuzzy.cpp src/refGeneParser.hpp src/polarUtility.hpp
	@echo "[Compile] perf_fuzzy"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) perf_testing/perf_fuzzy.cpp -o $(TPATH)perf_fuzzy.o

$(TPATH)refGeneParser.o : src/refGeneParser.hpp src/refGeneParser.cpp
	@echo "[Compile] RefGeneParser"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/refGeneParser.cpp -o $(TPATH)refGeneParser.o

$(TPATH)perf_bppFuzzy : $(TPATH)perf_bppFuzzy.o \
	$(TPATH)createTPset.o \
	$(TPATH)createTNset.o \
	$(TPATH)hgvsParser.o \
	lib/polarUtility.o \
	$(TPATH)readKnownPolyA.o \
	$(TPATH)refGeneParser.o \
	$(TPATH)utr3Finder.o \
	$(TPATH)utr3FinderFuzzy.o \
	$(TPATH)utr3FinderNaive.o \
	$(TPATH)seqStruct.o \
	$(TPATH)bppPredictFuzzy.o
	@echo "[Link] perf_bppFuzzy"
	@$(CC) $(BOOST) $(SEQAN) $^ $(LIBPATH) $(LFLAGS) $(LIBS) -o $(TPATH)perf_bppFuzzy

$(TPATH)perf_bppFuzzy.o : perf_testing/perf_bppFuzzy.hpp perf_testing/perf_bppFuzzy.cpp src/polarUtility.hpp
	@echo "[Compile] perf_bppFuzzy"
	@$(CC) $(BOOST) $(SEQAN) $(VIENNA) $(LIBPATH) $(CFLAGS) $(LIBS) perf_testing/perf_bppFuzzy.cpp -o $(TPATH)perf_bppFuzzy.o

.PHONY : clean
clean :
	@echo "[Delete] object and binary files"
	@rm $(TPATH)polar $(TPATH)perf_fuzzy $(TPATH)uTests $(TPATH)acc_test $(TPATH)*.o 
 

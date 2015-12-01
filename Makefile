CC 		= g++
CFLAGS		= -c -Wall -Wextra -pedantic -std=c++11 -O2 -g
LFLAGS		= -Wall
BOOST		= -I $(HOME)/boost_1_58_0 
SEQAN		= -I $(HOME)/seqan-library-2.0.0/include
LIBPATH		= -Llib/
CXX		= /usr/bin/gcc
LIBS		= -lboost_regex -lboost_filesystem -lboost_system -lboost_iostreams
LIBSTEST	= -lboost_unit_test_framework -lboost_filesystem -lboost_system -lboost_iostreams
TPATH		= bin/

OBJS		= $(TPATH)buildIndexFile.o \
	$(TPATH)fastaReader.o \
	$(TPATH)hgvsParser.o \
	$(TPATH)jannovarVcfParser.o \
	$(TPATH)knownGeneParser.o \
	$(TPATH)buildSeqStruct.o \
	$(TPATH)readTranscripts.o \
	$(TPATH)seqStruct.o \
	$(TPATH)utr3Finder.o \
	$(TPATH)utr3FinderFuzzy.o \
	$(TPATH)utr3FinderNaive.o 


.PHONY : all
all : $(TPATH)polar $(TPATH)uTests $(TPATH)acc_test $(TPATH)perf_fuzzy

$(TPATH)polar : $(TPATH)polar.o $(OBJS) 
	@echo "[Link] polar"
	@$(CC) $(BOOST) $(SEQAN) $^ $(LIBPATH) $(LFLAGS) $(LIBS) -o $(TPATH)polar 

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
	src/jannovarVcfParser.hpp \
	src/knownGeneParser.hpp \
	src/buildSeqStruct.hpp \
	src/readTranscripts.hpp \
	src/seqStruct.hpp \
	src/utr3Finder.hpp \
	src/utr3FinderFuzzy.hpp \
	src/utr3FinderNaive.hpp 
	@echo "[Compile] polar"
	@$(CC) $(BOOST) $(SEQAN) $(CFLAGS) src/polar.cpp -o $(TPATH)polar.o

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

#targets for the performance tests
$(TPATH)acc_test : $(TPATH)acc_test.o \
	$(TPATH)readKnownPolyA.o \
	$(TPATH)hgvsParser.o \
	$(TPATH)utr3Finder.o \
	$(TPATH)utr3FinderFuzzy.o \
	$(TPATH)utr3FinderNaive.o \
	$(TPATH)seqStruct.o
	@echo "[Link] acc_test"
	@$(CC) $(BOOST) $(SEQAN) $^ $(LIBPATH) $(LFLAGS) $(LIBS) -o $(TPATH)acc_test

$(TPATH)acc_test.o : perf_testing/acc_test.hpp perf_testing/acc_test.cpp perf_testing/readKnownPolyA.hpp
	@echo "[Compile] acc_test"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) perf_testing/acc_test.cpp -o $(TPATH)acc_test.o

$(TPATH)readKnownPolyA.o : perf_testing/readKnownPolyA.hpp perf_testing/readKnownPolyA.cpp 
	@echo "[Compile] readKnownPolyA"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) perf_testing/readKnownPolyA.cpp -o $(TPATH)readKnownPolyA.o

$(TPATH)perf_fuzzy : $(TPATH)perf_fuzzy.o \
	$(TPATH)readKnownPolyA.o \
	$(TPATH)hgvsParser.o \
	$(TPATH)refGeneParser.o \
	$(TPATH)utr3Finder.o \
	$(TPATH)utr3FinderFuzzy.o \
	$(TPATH)utr3FinderNaive.o \
	$(TPATH)seqStruct.o
	@echo "[Link] perf_fuzzy"
	@$(CC) $(BOOST) $(SEQAN) $^ $(LIBPATH) $(LFLAGS) $(LIBS) -o $(TPATH)perf_fuzzy

$(TPATH)perf_fuzzy.o : perf_testing/perf_fuzzy.hpp perf_testing/perf_fuzzy.cpp
	@echo "[Compile] perf_fuzzy"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) perf_testing/perf_fuzzy.cpp -o $(TPATH)perf_fuzzy.o

$(TPATH)refGeneParser.o : src/refGeneParser.hpp src/refGeneParser.cpp
	@echo "[Compile] RefGeneParser"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) src/refGeneParser.cpp -o $(TPATH)refGeneParser.o



.PHONY : clean
clean :
	@echo "[Delete] object and binary files"
	@rm $(TPATH)polar $(TPATH)uTests $(TPATH)acc_test $(TPATH)*.o 
 

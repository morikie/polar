CC 		= g++
CFLAGS		= -c -Wall -Wextra -pedantic -std=c++11 -O2 -g
LFLAGS		= -Wall
INCLUDE		= -I $(HOME)/boost_1_58_0 -I $(HOME)/seqan-library-2.0.0/include
LIBPATH		= -Llib/
CXX		= /usr/bin/gcc
LIBS		= -lboost_regex -lboost_filesystem -lboost_system
LIBSTEST	= -lboost_unit_test_framework -lboost_filesystem -lboost_system
TPATH		= bin/

OBJS		= $(TPATH)fastaReader.o \
	$(TPATH)utr3Finder.o \
	$(TPATH)seqStruct.o \
	$(TPATH)hgvsParser.o \
	$(TPATH)readTranscripts.o \
	$(TPATH)knownGeneParser.o \
	$(TPATH)jannovarVcfParser.o \
	$(TPATH)hgvsParser.o \
	$(TPATH)readSeqStruct.o 


.PHONY : all
all : $(TPATH)polar $(TPATH)uTests

$(TPATH)polar : $(TPATH)polar.o $(OBJS) 
	@echo "[Link] polar"
	@$(CC) $(INCLUDE) $^ $(LIBPATH) $(LFLAGS) $(LIBS) -o $(TPATH)polar 

$(TPATH)uTests : $(TPATH)uTests.o $(OBJS) 
	@echo "[Link] uTests"
	@$(CC) $^ $(LIBPATH) $(LFLAGS) $(LIBSTEST) -o bin/uTests

$(TPATH)uTests.o : unit_testing/uTests.cpp
	@echo "[Compile] uTests"
	@$(CC) $(INCLUDE) $(LIBPATH) $(CFLAGS) unit_testing/uTests.cpp -o $(TPATH)uTests.o

$(TPATH)polar.o : src/polar.hpp src/polar.cpp src/utr3Finder.hpp src/knownGeneParser.hpp src/jannovarVcfParser.hpp src/readTranscripts.hpp src/seqStruct.hpp
	@echo "[Compile] polar"
	@$(CC) $(INCLUDE) $(CFLAGS) src/polar.cpp -o $(TPATH)polar.o

$(TPATH)fastaReader.o : src/fastaReader.hpp src/fastaReader.cpp
	@echo "[Compile] FastaReader"
	@$(CC) $(INCLUDE) $(LIBPATH) $(CFLAGS) $(LIBS) src/fastaReader.cpp -o $(TPATH)fastaReader.o

$(TPATH)utr3Finder.o : src/utr3Finder.hpp src/utr3Finder.cpp src/hgvsParser.hpp src/seqStruct.hpp
	@echo "[Compile] UTR3MutationFinder"
	@$(CC) $(INCLUDE) $(LIBPATH) $(CFLAGS) $(LIBS) src/utr3Finder.cpp -o $(TPATH)utr3Finder.o

$(TPATH)seqStruct.o : src/seqStruct.hpp src/seqStruct.cpp src/hgvsParser.hpp
	@echo "[Compile] SeqStruct"
	@$(CC) $(INCLUDE) $(LIBPATH) $(CFLAGS) $(LIBS) src/seqStruct.cpp -o $(TPATH)seqStruct.o

$(TPATH)readTranscripts.o : src/readTranscripts.hpp src/readTranscripts.cpp
	@echo "[Compile] ReadTranscripts"
	@$(CC) $(INCLUDE) $(LIBPATH) $(CFLAGS) $(LIBS) src/readTranscripts.cpp -o $(TPATH)readTranscripts.o

$(TPATH)knownGeneParser.o : src/knownGeneParser.hpp src/knownGeneParser.cpp
	@echo "[Compile] KnownGeneParser"
	@$(CC) $(INCLUDE) $(LIBPATH) $(CFLAGS) $(LIBS) src/knownGeneParser.cpp -o $(TPATH)knownGeneParser.o

$(TPATH)jannovarVcfParser.o : src/jannovarVcfParser.hpp src/jannovarVcfParser.cpp
	@echo "[Compile] JannovarVcfParser"
	@$(CC) $(INCLUDE) $(LIBPATH) $(CFLAGS) $(LIBS) src/jannovarVcfParser.cpp -o $(TPATH)jannovarVcfParser.o

$(TPATH)hgvsParser.o : src/hgvsParser.hpp src/hgvsParser.cpp
	@echo "[Compile] HGVS Parser"
	@$(CC) $(INCLUDE) $(LIBPATH) $(CFLAGS) $(LIBS) src/hgvsParser.cpp -o $(TPATH)hgvsParser.o

$(TPATH)readSeqStruct.o : src/readSeqStruct.hpp src/readSeqStruct.cpp
	@echo "[Compile] readSeqStruct"
	@$(CC) $(INCLUDE) $(LIBPATH) $(CFLAGS) $(LIBS) src/readSeqStruct.cpp -o $(TPATH)readSeqStruct.o


.PHONY : clean
clean :
	@echo "[Delete] object and binary files"
	@rm $(TPATH)polar $(TPATH)uTests $(TPATH)*.o
 

CC 		= g++
CFLAGS		= -Wall -Wextra -g -pedantic -std=c++11
CINCLUDE	= -I /home/morikie/boost_1_58_0 -I /home/morikie/seqan-library-2.0.0
CLIBPATH	= -Llib/
CXX		= /usr/bin/gcc
CLIBS		= -lboost_regex -lboost_filesystem -lboost_system
CLIBSTEST	= -lboost_unit_test_framework -lboost_filesystem -lboost_system
TPATH		= bin/

 
.PHONY : all
all : $(TPATH)polar $(TPATH)uTests

$(TPATH)polar : $(TPATH)polar.o $(TPATH)fastaReader.o $(TPATH)mutationFinder.o
	@echo "[Link] polar"
	@$(CC) $(CINCLUDE) -o $(TPATH)polar $^ $(CLIBPATH) $(CLIBS) $(CFLAGS)

$(TPATH)uTests : $(TPATH)uTests.o $(TPATH)fastaReader.o $(TPATH)mutationFinder.o
	@echo "[Link] uTests"
	@$(CC) -o bin/uTests $^ $(CLIBPATH) $(CLIBSTEST) $(CFLAGS) 

$(TPATH)uTests.o : unit_testing/uTests.cpp
	@echo "[Compile] uTests"
	@$(CC) $(CINCLUDE) $(CLIBPATH) -c unit_testing/uTests.cpp -o $(TPATH)uTests.o $(CFLAGS)

$(TPATH)polar.o : src/polar.cpp src/polar.hpp
	@echo "[Compile] polar"
	@$(CC) -c $< -o $(TPATH)polar.o $(CFLAGS)

$(TPATH)fastaReader.o : src/fastaReader.cpp src/fastaReader.hpp
	@echo "[Compile] FastaReader"
	@$(CC) $(CINCLUDE) -c $< -o $(TPATH)fastaReader.o $(CLIBPATH) $(CLIBS) $(CFLAGS)

$(TPATH)mutationFinder.o : src/mutationFinder.cpp src/mutationFinder.hpp
	@echo "[Compile] MutationFinder"
	$(CC) $(CINCLUDE) -c $< -o $(TPATH)mutationFinder.o $(CLIBPATH) $(CLIBS) $(CFLAGS)

.PHONY : clean
clean :
	@echo "[Delete] object and binary files"  
	@rm $(TPATH)polar $(TPATH)uTests $(TPATH)*.o
 

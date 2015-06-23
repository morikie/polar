CC 		= g++
CFLAGS		= -Wall -g -pedantic
CINCLUDE	= -I /home/morikie/boost_1_58_0
CLIBPATH	= -Llib/
CXX		= /usr/bin/gcc
CLIBS		= -lboost_regex -lboost_filesystem -lboost_system
CLIBSTEST	= -lboost_unit_test_framework -lboost_filesystem -lboost_system
TPATH		= bin/

 
.PHONY : all
all : $(TPATH)polar $(TPATH)uTests

$(TPATH)polar : $(TPATH)polar.o $(TPATH)fastaReader.o
	@echo "[Link] polar"
	$(CC) $(CINCLUDE) -o $(TPATH)polar $^ $(CLIBPATH) $(CLIBS) $(CFLAGS)

$(TPATH)uTests : $(TPATH)uTests.o $(TPATH)fastaReader.o
	@echo "[Link] uTests"
	$(CC) -o bin/uTests $^ $(CLIBPATH) $(CLIBSTEST) 

$(TPATH)uTests.o : unit_testing/uTests.cpp
	@echo "[Compile] uTests"
	$(CC) $(CINCLUDE) $(CLIBPATH) -c unit_testing/uTests.cpp -o $(TPATH)uTests.o $(CFLAGS)

$(TPATH)polar.o : src/polar.cpp src/polar.hpp
	@echo "[Compile] polar"
	$(CC) -c src/polar.cpp -o $(TPATH)polar.o

$(TPATH)fastaReader.o : src/fastaReader.cpp src/fastaReader.hpp
	@echo "[Compile] FastaReader"
	@echo "$<"
	$(CC) $(CINCLUDE) -c $< -o $(TPATH)fastaReader.o $(CLIBPATH) $(CLIBS)

.PHONY : clean
clean : 
	rm $(TPATH)polar $(TPATH)uTests $(TPATH)*.o
 

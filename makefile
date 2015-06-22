CC 		= g++
CFLAGS		= -Wall -g -pedantic 
CINCLUDE	= -I /home/morikie/boost_1_58
CLIBPATH	= -Llib/
CXX		= /usr/bin/gcc
CLIBS		= -lboost_regex
TPATH		= bin/

 
.PHONY : all
all : $(TPATH)polar $(TPATH)uTests

$(TPATH)polar : $(TPATH)polar.o
	$(CC) -o $(TPATH)polar $(TPATH)polar.o $(CLIBPATH) $(CLIBS) $(CFLAGS) $(CINCLUDE)

$(TPATH)uTests : $(TPATH)uTests.o
	$(CC) -o bin/uTests $(TPATH)uTests.o 

$(TPATH)uTests.o : unit_testing/uTests.cpp
	$(CC) -c unit_testing/uTests.cpp -o $(TPATH)uTests.o $(CFLAGS) $(CINCLUDE)

$(TPATH)polar.o : src/polar.cpp src/polar.hpp
	$(CC) -c src/polar.cpp -o bin/polar.o

.PHONY : clean
clean : 
	rm bin/*
 

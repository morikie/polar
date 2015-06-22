#define BOOST_TEST_DYN_LINK
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <boost/bind.hpp>
#include <boost/test/included/unit_test.hpp>
#include "../src/hello.hpp"

using namespace boost::unit_test;
using namespace std;

void testFunction(int i, int j) {
	BOOST_CHECK( i == j );
}


bool initFunction() {
	Auto golf;
	golf.hp = 20;

	framework::master_test_suite().add(BOOST_TEST_CASE( boost::bind(&testFunction, golf.getPs(), 50)));

	return true;
}

int main (int argc, char * argv[]) {
	
	return ::boost::unit_test::unit_test_main(&initFunction, argc, argv);
}


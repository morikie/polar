#include <cstdlib>
#include <iostream>
#include <vector>
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>
#include "polar.hpp"

using namespace std;


const set<string> utr3RecogHexamers = boost::assign::list_of
	("AATAAA")
	("ATTAAA")
	("TATAAA")
	("AGTAAA")
	("AAGAAA")
	("AATATA")
	("AATACA")
	("CATAAA")
	("GATAAA")
	("AATGAA")
	("TTTAAA")
	("ACTAAA")
	("AATAGA")
;


int main (int args, char * argv[]) {
	int numbers[] = {1,2,3,4,5};
	vector<int> test (numbers, numbers + 5);
	BOOST_FOREACH (int i, test) {
		cout << i << " ";
	}
	cout << endl;
	cout << "Hello World" << endl;
	boost::regex e("mau.");
	string str1 = "marder";
	string str2 = "maus";
	cout << regex_match(str1, e) << endl;
	cout << regex_match(str2, e) << endl;
	return EXIT_SUCCESS;
}

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <set>
#include <iterator>
#include <algorithm>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/assign.hpp>
#include <cmath>
#include <chrono>

using std::string;
using std::ostream;
using std::vector;
using std::set;
using std::map;
using std::tuple;
using std::ostringstream;
using std::to_string;
using std::istream;
using std::ostream;
using std::cin;
using std::get;
using std::istringstream;
using std::stoi;
using std::stod;
using std::min;

template<typename T, typename C>
ostream& operator<<(ostream& os, const set<T, C>& c) {
	for (auto const &i : c)
		os << i << " ";
	return os;
}

template<typename T>
ostream& operator<<(ostream& os, const vector<T>& c) {
	for (auto const &i : c)
		os << i << " ";
	return os;
}

template<typename T, typename K>
ostream& operator<<(ostream& os, const map<T, K>& c) {
	for (auto &i : c)
		os << i.first << ":" << i.second << " ";
	return os;
}

template<typename T>
string array2string(const T* a, const T* b) {
	ostringstream os;
	for (const T* i=a; i<b; i++)
		os << *i << " ";
	return os.str();
}

#define LOGGER(X)
#ifdef ENABLE_LOGGER
#undef LOGGER
#define LOGGER(X) X
#define logger std::cerr
#endif
using std::endl;


long long power(int a, int b, int mod) {
	if (b == 0)
		return 1;
	long long p = power(a, b/2, mod);
	p = (p * p) % mod;
	if (b % 2 == 1)
		p = p * a;
	return p % mod;
}




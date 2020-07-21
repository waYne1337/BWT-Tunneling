/*
 * tfm_index_invert.cpp for Edge minimization in de Bruijn graphs
 * Copyright (c) 2019 Uwe Baier, Pascal Weber All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <iostream>
#include <map>
#include <deque>
#include <utility>
#include <vector>

#include "tfm_index.hpp"

using namespace std;
using namespace sdsl;

typedef typename sdsl::int_vector<>::size_type size_type;

void printUsage( char **argv ) {
	cerr << "USAGE: " << argv[0] << " TFMFILE" << endl;
	cerr << "TFMFILE:" << endl;
	cerr << "  File where to store the serialized trie" << endl;
};

int main(int argc, char **argv) {
	//check parameters
	if (argc != 2) {
		printUsage( argv );
		cerr << "At least 1 parameter expected" << endl;
		return 1;
	}
	
	//load tunneled fm index
	tfm_index<> tfm;
	load_from_file( tfm, argv[1] );

	//reconstruct original string using tunneled fm index
	char *S = new char[tfm.size()];
	S[tfm.size()-1] = '\0';

	auto p = tfm.end();
	for (size_type i = 1; i < tfm.size(); i++) {
		S[tfm.size() - i - 1] = (char)tfm.backwardstep( p );
	}

	//print result
	cout << S;
	return 0;
}

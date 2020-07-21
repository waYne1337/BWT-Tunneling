/*
 * dbg_edgespectrum.cpp for BWT Tunneling
 * Copyright (c) 2020 Uwe Baier All Rights Reserved.
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
//ORIGINALLY COMES FROM:
/*
 * dbg_edgespectrum.cpp for Edge minimization in de Bruijn graphs
 * Copyright (c) 2019 Uwe Baier, Pascal Weber All Rights Reserved.
 */

#include <iostream>
#include <stdlib.h>
#include <string>

#include <sdsl/int_vector.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/util.hpp>

#include "dbg_algorithms.hpp"

using namespace std;
using namespace sdsl;

typedef typename sdsl::int_vector<>::size_type size_type;

void printUsage( const char *command ) {
	cerr << command << " k file" << endl;
	cerr << "\tk: maximal order s.t. amount of edges from edge-reduced DBGS of order [1,k] are listed" << endl;
	cerr << "\tfile: input file using byte alphabet (mustn't contain nullbytes)" << endl;
}

int main(int argc, char **argv ) {
	if (argc != 3) {
		printUsage(argv[0]);
		return 1;
	}
	//fetch args
	string infile = argv[2];
	size_type k = atoi( argv[1] );
	csa_wt<wt_blcd<>,0xFFFFFFFF,0xFFFFFFFF> csa;
	cache_config config(true, "./", util::basename(infile) );
	construct( csa, argv[2], config, 1 );
	//check correctness
	if (k == 0 || k > csa.size()) {
		printUsage( argv[0] );
		cerr << "k must be between 1 and size of file" << endl;
		return 1;
	}

	//compute edge spectrum
	auto ES = dbg_algorithms::dbg_edgespectrum( csa, k );

	//output results
	for (size_type i = 0; i < k; i++) {
		cout << "k" << (i+1) << "\t" << ES[i] << endl;
	}
	return 0;
}
	
	
	

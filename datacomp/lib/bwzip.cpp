/*
 * bwzip.cpp for BWT Tunneling
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
 * ui.cpp for bwt tunneling
 * Copyright (c) 2017 Uwe Baier All Rights Reserved.
 */

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string.h>

#include "bwt_compressor.hpp"

//tunnel planning strategies
#include "tp_strategy_none.hpp"
#include "tp_strategy_hirsch.hpp"
#include "tp_strategy_greedy.hpp"
#include "tp_strategy_greedy_update.hpp"
#include "tp_strategy_bestp.hpp"

//post stages
#include "bcm_poststage.hpp"
#include "bw94_poststage.hpp"

using namespace std;

enum bwt_tunnel_strategy{ NONE, HIRSCH, GREEDY, GREEDY_UPDATE, BESTP };
enum bwt_post_stage{ BW94, BCM };

//forward declarations
template<typename t_tunnel_strat>
int bw_compress( istream &in, ostream &out, bwt_post_stage post_stage, bool informative );
template<typename t_tunnel_strat, typename t_post_stage>
int bw_compress( istream &in, ostream &out, bool informative );

template<typename t_tunnel_strat>
int bw_decompress( istream &in, ostream &out, bool informative );
template<typename t_tunnel_strat, typename t_post_stage>
int bw_decompress( istream &in, ostream &out, bool informative );

void printUsage( char **argv ) {
	cerr << "USAGE: " << argv[0] << " [OPTIONS] INFILE OUTFILE" << endl;
	cerr << "OPTIONS:" << endl;
	cerr << "  -d\tdecompress data (compression is default)." << endl;
	cerr << "    \tIf enabled, ignores all except of the -i options." << endl;
	cerr << "  -i\tEnable informative mode, printing additional information" << endl;
	cerr << "  -tstrat [STRATEGY]\ttunneling strategy to be used. Must be one of the following:" << endl;
	cerr << "                    \tnone : enable no tunneling" << endl;
	cerr << "                    \thirsch : hirsch tunnel planning strategy (default)" << endl;
	cerr << "                    \tgreedy : greedy strategy considering no side effects" << endl;
	cerr << "                    \tgreedy-update : greedy strategy considering negative side effects" << endl;
	cerr << "                    \tbest[PERCENT] : tunnel the best PERCENT of prefix intervals," << endl;
	cerr << "                    \t                e.g. best10 for the best 10 %" << endl;
	cerr << "  -pstage [PSTAGE]\tpost stages used for the compression of the BWT. Must be one of" << endl;
	cerr << "                  \tbw94 : compression scheme from 1994 using move-to-front transform," << endl;
	cerr << "                  \t       run-length encoding and source encoding, as described" << endl;
	cerr << "                  \t       by Mike Burrows and David J. Wheeler" << endl;
	cerr << "                  \tbcm :  compression using a bwt-optimized context mixer (default)" << endl;
	cerr << "                  \t       by Ilya Muravyov" << endl;
	cerr << "INFILE:" << endl;
	cerr << "  File to be compressed or decompressed if -d is set" << endl;
	cerr << "OUTFILE:" << endl;
	cerr << "  File to store the compressed (or uncompressed) data, see -d flag." << endl;
};

int main( int argc, char **argv ) {
	//analyse args
	bool compress = true; //compress or decompress
	bool informative = false; //informative mode
	string infile;
	string outfile;

	//check parameters
	if (argc < 3) {
		printUsage( argv );
		cerr << "At least 2 parameters expected" << endl;
		return 1;
	}

	//tunnel strategy
	bwt_tunnel_strategy tunnel_strategy = HIRSCH;

	//post bwt stage
	bwt_post_stage post_stage = BCM;

	//last option in arguments
	enum {NO, COMP, INF, TSTRAT, PSTAGE} last_option;
	last_option = NO;

	for (int i = 1; i < argc - 2; i++) { //analyze options
		switch (last_option) {
		case NO: //last options that require no additional parameter
		case COMP:
		case INF:
			if (strcmp(argv[i], "-d") == 0) { //decompress
				last_option = COMP;
				compress = false;
			}
			else if (strcmp(argv[i], "-i") == 0) {
				last_option = INF;
				informative = true;
			}
			else if (strcmp(argv[i], "-tstrat") == 0) {
				last_option = TSTRAT;
			}
			else if (strcmp(argv[i], "-pstage") == 0) {
				last_option = PSTAGE;
			}
			else {
				printUsage(argv);
				cerr << "Unknown option " << argv[i] << endl;
				return 1;
			}
			break;
		case TSTRAT: //determine tunneling strategy
			if (strcmp(argv[i], "none") == 0) {
				 tunnel_strategy = NONE;
			}
			else if (strcmp(argv[i], "hirsch") == 0) {
				 tunnel_strategy = HIRSCH;
			}
			else if (strcmp(argv[i], "greedy") == 0) {
				 tunnel_strategy = GREEDY;
			}
			else if (strcmp(argv[i], "greedy-update") == 0) {
				 tunnel_strategy = GREEDY_UPDATE;
			}
			else if (strncmp(argv[i], "best", 4) == 0) {
				int p = atoi(argv[i] + 4);
				if (p < 0 || p > 100) {
					cerr << "value for the best p percentage must be between 0 and 100" << endl;
					return 1;
				}
				tp_strategy_bestp::percentage = p;
				tunnel_strategy = BESTP;
			}
			else {
				printUsage(argv);
				cerr << "Unknown tunneling strategy option " << argv[i] << endl;
				return 1;
			}
			last_option = NO;
			break;
		case PSTAGE: //determine post BWT stage
			if (strcmp(argv[i], "bcm") == 0) {
				post_stage = BCM;
			}
			else if (strcmp(argv[i], "bw94") == 0) {
				post_stage = BW94;
			}
			else {
				printUsage(argv);
				cerr << "Unknown post stage option " << argv[i] << endl;
				return 1;
			}
			last_option = NO;
			break;
		}
	}
	infile = argv[argc-2];
	outfile = argv[argc-1];

	//open streams for infile and outfile
	ifstream fin{ infile };
	ofstream fout{ outfile, ofstream::out | ofstream::trunc };
	if (!fin) {
		printUsage(argv);
		cerr << "unable to open file \"" << infile << "\"" << endl;
		return 1;
	} else if (!fout) {
		printUsage(argv);
		cerr << "unable to open file \"" << outfile << "\"" << endl;
		return 1;
	}

	if (compress) {
		switch (tunnel_strategy) {
		case NONE:
			fout << "non" << endl;
			return bw_compress<tp_strategy_none>( fin, fout, post_stage, informative );
		case HIRSCH:
			fout << "hir" << endl;
			return bw_compress<tp_strategy_hirsch>( fin, fout, post_stage, informative );
		case GREEDY:
			fout << "grd" << endl;
			return bw_compress<tp_strategy_greedy>( fin, fout, post_stage, informative );
		case GREEDY_UPDATE:
			fout << "gdu" << endl;
			return bw_compress<tp_strategy_greedy_update>( fin, fout, post_stage, informative );
		case BESTP:
			fout << "bep" << endl;
			return bw_compress<tp_strategy_bestp>( fin, fout, post_stage, informative );
		}
	} else {
		//read first line of input and decide what to do
		string tstrat;
		getline( fin, tstrat );
		if (tstrat == "non") {
			return bw_decompress<tp_strategy_none>( fin, fout, informative );
		}
		else if (tstrat == "hir") {
			return bw_decompress<tp_strategy_hirsch>( fin, fout, informative );
		}
		else if (tstrat == "grd") {
			return bw_decompress<tp_strategy_greedy>( fin, fout, informative );
		}
		else if (tstrat == "gdu") {
			return bw_decompress<tp_strategy_greedy_update>( fin, fout, informative );
		}
		else if (tstrat == "bep") {
			return bw_decompress<tp_strategy_bestp>( fin, fout, informative );
		}
		printUsage( argv );
		cerr << "Unknown tunneling strategy " << tstrat << "in the encoding of " << infile << ", unable to decompress" << endl;
		return 1;
	}
	return 0;	
}

template<typename t_tunnel_strat>
int bw_compress( istream &in, ostream &out, bwt_post_stage post_stage, bool informative ) {
	//output post stage identifier
	switch (post_stage) {
	case BW94:
		out << "w94" << endl;
		return bw_compress<t_tunnel_strat,bw94_poststage>( in, out, informative );
	case BCM:
		out << "bcm" << endl;
		return bw_compress<t_tunnel_strat,bcm_poststage>( in, out, informative );
	}
	return 1;
}

template<typename t_tunnel_strat>
int bw_decompress( istream &in, ostream &out, bool informative ) {
	//read post stage
	string post_stage;
	getline( in, post_stage );
	if (post_stage == "w94") {
		return bw_decompress<t_tunnel_strat,bw94_poststage>( in, out, informative );
	}
	else if (post_stage == "bcm") {
		return bw_decompress<t_tunnel_strat,bcm_poststage>( in, out, informative );
	}
	else {
		cerr << "Unknown post stage " << post_stage << " used to compress file, unable to decompress" << endl;
		return 1;
	}
}

template<typename t_tunnel_strat, typename t_post_stage>
int bw_compress( istream &in, ostream &out, bool informative ) {
	bwt_compressor<t_tunnel_strat,t_post_stage> compressor;
	compressor.set_quiet( !informative );
	compressor.compress( in, out );
	return 0;
}

template<typename t_tunnel_strat, typename t_post_stage>
int bw_decompress( istream &in, ostream &out, bool informative ) {
	bwt_compressor<t_tunnel_strat,t_post_stage> compressor;
	compressor.set_quiet( !informative );
	compressor.decompress( in, out );
	return 0;
}

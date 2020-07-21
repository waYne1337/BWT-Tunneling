/*
 * trie_construct.cpp for BWT Tunneling
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

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <sdsl/config.hpp>
#include <sdsl/construct_config.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

#include "tries.hpp"

using namespace std;
using namespace sdsl;

typedef typename sdsl::int_vector<>::size_type size_type;

//construction dependent on trie type
template<class trie_t>
int construct_trie( const std::string &infile, const std::string &outfile, bool informative );

void printUsage( char **argv ) {
	cerr << "USAGE: " << argv[0] << " [OPTIONS] INFILE TRIEOUTFILE" << endl;
	cerr << "OPTIONS:" << endl;
	cerr << "  -i\tEnable informative mode, printing memory peaks (in bytes) and" << endl;
	cerr << "    \tconstruction timings (in milliseconds) during construction" << endl;
	cerr << "  -sa\tChoose suffix array construction algorithm. Must be followed by one of:" << endl;
	cerr << "     \tDIVSUFSORT      use divsufsort (fast but memory-intensive)" << endl;
	cerr << "     \tSE_SAIS         use a semi-external algorithm (slower but lower mem peak)" << endl;
	cerr << "  -ta\tChoose trie construction algorithm. Must be followd by one of:" << endl;
	cerr << "     \tXBWT        fast but memory-intensive" << endl;
	cerr << "     \tXBWT_SC     fast and low mem peak" << endl;
	cerr << "     \tXBWT_LW     slow and low mem peak" << endl;
	cerr << "     \tXBWT_LW_SC  slow and lowest mem peak" << endl;
	cerr << "     \tTXBWT       medium but memory-intensive, trie compression" << endl;
	cerr << "     \tTXBWT_SC    medium and low mem peak, trie compression" << endl;
	cerr << "INFILE:" << endl;
	cerr << "  File to construct trie from. strings must be seperated by newlines, nullbytes are permitted" << endl;
	cerr << "TRIEOUTFILE:" << endl;
	cerr << "  File where to store the serialized trie" << endl;
};

//little hack to get extra information from memory managemant
const format_type leet_format = (format_type)1337;
namespace sdsl {
	template<>
	void write_mem_log<leet_format>(ostream& out,const memory_monitor& m) {

		//get all memory events
		auto events = m.completed_events;
		std::sort( events.begin(), events.end() );
		auto e = events.begin();

		//scan events and detect peak and time of both suffix array and xbwt construction

		//scan for fm index construction
		int64_t fm_mem_peak = 0;
		auto fm_start = m.start_log;
		auto fm_end = fm_start;		
		while (e != events.end() && e->name != "TRIECONSTRUCT") {
			for (auto alloc : e->allocations) {
				fm_mem_peak = std::max(fm_mem_peak, alloc.usage);
				fm_end = alloc.timestamp;
			}
			++e;
		}

		//scan for trie construction
		int64_t trie_mem_peak = 0;
		auto trie_start = fm_end;
		auto trie_end = trie_start;
		while (e != events.end()) {
			for (auto alloc : e->allocations) {
				trie_mem_peak = std::max(trie_mem_peak, alloc.usage);
				trie_end = alloc.timestamp;
			}
			++e;
		}

		//print results
		out << "fm_mem_peak\t" << fm_mem_peak << endl;
		out << "fm_time\t" << chrono::duration_cast<chrono::milliseconds>(fm_end-fm_start).count() << endl;

		out << "trie_mem_peak\t" << trie_mem_peak << endl;
		out << "trie_time\t" << chrono::duration_cast<chrono::milliseconds>(trie_end-trie_start).count() << endl;
	};
};

//main program
int main(int argc, char **argv) {
	//set default configuration
	construct_config::byte_algo_sa = LIBDIVSUFSORT;
	trie_construct_config::algo = XBWT_SC;
	bool informative = false; //informative mode
	string infile = "Makefile";
	string outfile = "Makefile.xbwt";

	//check parameters
	if (argc < 3) {
		printUsage( argv );
		cerr << "At least 2 parameters expected" << endl;
		return 1;
	}
	enum { SA, TA, IN, NO } last_option; //enumeration for last option
	last_option = NO;
	for (int i = 1; i < argc - 2; i++) { //analyze options
		switch (last_option) {
		case NO:
		case IN: //last option does not require a parameter, read next option
			if (strcmp(argv[i], "-i") == 0) { //enable informative mode
				last_option = IN;
				informative = true;
			}
			else if (strcmp(argv[i], "-sa") == 0) {
				last_option = SA;
			}
			else if (strcmp(argv[i], "-ta") == 0) {
				last_option = TA;
			}
			else {
				printUsage(argv);
				cerr << "Unknown option " << argv[i] << endl;
				return 1;
			}
			break;
		case SA: //choose suffix array construction algorithm
			if (strcmp(argv[i], "DIVSUFSORT") == 0) {
				construct_config::byte_algo_sa = LIBDIVSUFSORT;
			}
			else if (strcmp(argv[i], "SE_SAIS") == 0) {
				construct_config::byte_algo_sa = SE_SAIS;
			}
			else {
				printUsage( argv );
				cerr << "Unknown suffix array construction algorithm " << argv[i] << endl;
				return 1;
			}
			last_option = NO;
			break;
		case TA: //choose trie construction algorithm
			if (strcmp(argv[i], "XBWT") == 0) {
				trie_construct_config::algo = XBWT;
			}
			else if (strcmp(argv[i], "XBWT_SC") == 0) {
				trie_construct_config::algo = XBWT_SC;
			}
			else if (strcmp(argv[i], "XBWT_LW") == 0) {
				trie_construct_config::algo = XBWT_LW;
			}
			else if (strcmp(argv[i], "XBWT_LW_SC") == 0) {
				trie_construct_config::algo = XBWT_LW_SC;
			}
			else if (strcmp(argv[i], "TXBWT") == 0) {
				trie_construct_config::algo = TXBWT;
			}
			else if (strcmp(argv[i], "TXBWT_SC") == 0) {
				trie_construct_config::algo = TXBWT_SC;
			}
			else {
				printUsage( argv );
				cerr << "Unknown trie construction algorithm " << argv[i] << endl;
				return 1;
			}
			last_option = NO;
			break;
		}
	}
	infile = argv[argc-2];
	outfile = argv[argc-1];

	switch (trie_construct_config::algo) {
	case XBWT:
	case XBWT_SC:
	case XBWT_LW:
	case XBWT_LW_SC:
		return construct_trie<trie_xbwt<>>( infile, outfile, informative );
	case TXBWT:
	case TXBWT_SC:
		return construct_trie<trie_txbwt<>>( infile, outfile, informative );
	}
}

template<class trie_t>
int construct_trie( const std::string &infile, const std::string &outfile, bool informative ) {
	//set up a config
	cache_config config(true, "./", util::basename(infile) );
	trie_t trie;

	//and run construction
	memory_monitor::start();
	construct( trie, infile, config, 1 );
	memory_monitor::stop();

	//save trie
	osfstream out(outfile, std::ios::binary | std::ios::trunc | std::ios::out);
	out << trie_trait<trie_t>::name << endl; //add typeinfo
	if (!out) {
		cerr << "Unable to store trie to file " << outfile << endl;
		return 1;
	}
	auto num_bytes = serialize(trie,out);
	out.close();

	if (informative) { //output construction information
		//print additional information
		cout << "strings_length\t" << trie.strings_length() << endl;
		cout << "num_strings\t" << trie.num_leaves() << endl;
		cout << "trie_size\t" << num_bytes << endl;

		memory_monitor::write_memory_log<leet_format>(cout);
	}
	return 0;
}

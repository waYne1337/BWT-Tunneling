/*
 * mp_search.cpp for BWT Tunneling
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
//ORIGNINALLY COMES FROM:
/*
 * mpsearch.cpp for trickier xbwt tricks
 * Copyright (c) 2018 Stefan Stauß, Uwe Baier All Rights Reserved.
 */
#include <algorithm>
#include <iostream>
#include <string>
#include <string.h>

#include <sdsl/sfstream.hpp>
#include <sdsl/memory_management.hpp>

using namespace std;

#include "tries.hpp"

using namespace std;
using namespace sdsl;

template<class trie_t>
void mp_search( isfstream &trie_in, bool informative );

void printUsage( char **argv ) {
	cerr << "USAGE: " << argv[0] << " [OPTIONS] TRIEFILE" << endl;
	cerr << "DESCRIPTION: Program searches for all occurrences of all listed patterns using stdin" << endl;
	cerr << "  Output is written to stdout and consists of a list of position pairs," << endl;
	cerr << "  where the first entry is the ending position of the found pattern, " << endl;
	cerr << "  and the second entry is the line number of the pattern in " << endl;
	cerr << "  the file used for trie construction." << endl;
	cerr << "OPTIONS:" << endl;
	cerr << "  -i\tEnable informative mode, printing memory peak (in bytes) and" << endl;
	cerr << "    \tsearch timing (in milliseconds) during search" << endl;
	cerr << "TRIEFILE: a file containing a serialized trie used for multi-pattern search" << endl;
	cerr << "NOTE: to work correct, trie should not contain strings which are proper substrings" << endl;
	cerr << "  of other contained strings, use `create_trie_input.x` to ensure this" << endl;
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

		//scan for timing
		int64_t mem_peak = 0;
		auto start = m.start_log;
		auto end = start;
		while (e != events.end()) {
			for (auto alloc : e->allocations) {
				mem_peak = std::max(mem_peak, alloc.usage);
				end = alloc.timestamp;
			}
			++e;
		}

		//print results
		out << "mpsearch_mem_peak\t" << mem_peak << endl;
		out << "mpsearch_time\t" << chrono::duration_cast<chrono::milliseconds>(end-start).count() << endl;
	};
};

//main program

int main( int argc, char **argv ) {
	//check arguments
	bool informative = false;
	string triefile;

	switch (argc) {
	case 2:
		triefile = argv[1];
		break;
	case 3:
		if (argv[1] != string("-i")) {
			printUsage( argv );
			return 1;
		}
		informative = true;
		triefile = argv[2];
		break;
	default:
		printUsage( argv );
		return 1;
	}
		
	//load serialized trie
	isfstream in(triefile, std::ios::binary | std::ios::in);
	if (!in) {
		cerr << "Unable to open file " << triefile << endl;
		return 1;
	}
	//check type of trie and run search
	string trietype;
	getline( in, trietype );
	if (trietype == trie_trait<trie_xbwt<>>::name) {
		mp_search<trie_xbwt<>>( in, informative );
	} else if (trietype == trie_trait<trie_txbwt<>>::name) {
		mp_search<trie_txbwt<>>( in, informative );
	} else {
		printUsage( argv );
		cerr << "Unknown trie format in file " << triefile << endl;
		return 1;
	}
	return 0;
}

template<class trie_t>
void mp_search( isfstream &trie_in, bool informative ) {
	trie_t trie;

	size_t i = -1; //current input position
	size_t em = 0; //number of exact matches
	memory_monitor::start();
	{
		auto event = memory_monitor::event("AHOCORASICK");
		trie.load( trie_in );

		//run aho corasick algorithm
		unsigned char c = '\0';
		auto v = trie.root();
		auto edges_v = trie.get_empty_edge_range();
		while (cin) {
			//navigate to next node
			auto edge = trie.find_edge( edges_v, c );
			if (edge.first != '\0') { //edge to an inner node found
				v = trie.follow_edge( edge, v );
				i++;
				c = cin.get();
			} else {
				if (v == trie.root()) {
					i++;
					c = cin.get();
				} else {
					v = trie.failure_link( v );
				}
			}
			edges_v = trie.get_edge_range( v );

			//check if a record is found
			size_t record = trie.has_leaf( edges_v );
			if (record != trie.num_leaves()) {
				cout << i << "\t" << record << endl;
				++em;
			}
		}
	}
	memory_monitor::stop();

	if (informative) { //output further information
		cout << "input_size\t" << i << endl;
		cout << "trie_strings_length\t" << trie.strings_length() << endl;
		cout << "trie_num_strings\t" << trie.num_leaves() << endl;
		cout << "trie_size\t" << size_in_bytes( trie ) << endl;
		cout << "num_exact_matches\t" << em << endl;

		memory_monitor::write_memory_log<leet_format>(cout);
	}
}

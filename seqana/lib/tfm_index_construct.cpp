/*
 * tfm_index_construct.cpp for BWT Tunneling
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
 * tfm_index_construct.cpp for Edge minimization in de Bruijn graphs
 * Copyright (c) 2019 Uwe Baier, Pascal Weber All Rights Reserved.
 */

#include <iostream>
#include <map>
#include <deque>
#include <utility>
#include <vector>

#include <sdsl/util.hpp>

#include "tfm_index.hpp"

using namespace std;
using namespace sdsl;

void printUsage( char **argv ) {
	cerr << "USAGE: " << argv[0] << " [OPTIONS] INFILE TFMOUTFILE" << endl;
	cerr << "OPTIONS:" << endl;
	cerr << "  -i\tEnable informative mode, printing memory peaks (in bytes) and" << endl;
	cerr << "    \tconstruction timings (in milliseconds) during construction" << endl;
	cerr << "  -sa\tChoose suffix array construction algorithm. Must be followed by one of:" << endl;
	cerr << "     \tDIVSUFSORT   use divsufsort (fast but memory-intensive)" << endl;
	cerr << "     \tSE_SAIS         use a semi-external algorithm (slower but lower mem peak)" << endl;
	cerr << "INFILE:" << endl;
	cerr << "  File to construct tunneled FM index from, nullbytes are permitted" << endl;
	cerr << "TFMOUTFILE:" << endl;
	cerr << "  File where to store the serialized tunneled fm index" << endl;
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
		while (e != events.end() && e->name != "FINDMINDBG") {
			for (auto alloc : e->allocations) {
				fm_mem_peak = std::max(fm_mem_peak, alloc.usage);
				fm_end = alloc.timestamp;
			}
			++e;
		}

		//scan for min dbg algorithm
		int64_t dbg_mem_peak = 0;
		auto dbg_start = fm_end;
		auto dbg_end = dbg_start;		
		while (e != events.end() && e->name != "TFMINDEXCONSTRUCT") {
			for (auto alloc : e->allocations) {
				dbg_mem_peak = std::max(dbg_mem_peak, alloc.usage);
				dbg_end = alloc.timestamp;
			}
			++e;
		}

		//scan for tunneled fm index construction
		int64_t tfm_mem_peak = 0;
		auto tfm_start = dbg_end;
		auto tfm_end = tfm_start;
		while (e != events.end()) {
			for (auto alloc : e->allocations) {
				tfm_mem_peak = std::max(tfm_mem_peak, alloc.usage);
				tfm_end = alloc.timestamp;
			}
			++e;
		}

		//print results
		out << "fm_mem_peak\t" << fm_mem_peak << endl;
		out << "fm_time\t" << chrono::duration_cast<chrono::milliseconds>(fm_end-fm_start).count() << endl;

		out << "min_dbg_mem_peak\t" << dbg_mem_peak << endl;
		out << "min_dbg_time\t" << chrono::duration_cast<chrono::milliseconds>(dbg_end-dbg_start).count() << endl;

		out << "tfm_mem_peak\t" << tfm_mem_peak << endl;
		out << "tfm_time\t" << chrono::duration_cast<chrono::milliseconds>(tfm_end-tfm_start).count() << endl;
	};
};

int main(int argc, char **argv) {
	//set default configuration
	construct_config::byte_algo_sa = LIBDIVSUFSORT;
	bool informative = false; //informative mode
	string infile = "Makefile";
	string outfile = "Makefile.tfm";

	//check parameters
	if (argc < 3) {
		printUsage( argv );
		cerr << "At least 2 parameters expected" << endl;
		return 1;
	}
	enum { SA, IN, NO } last_option; //enumeration for last option
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
		}
	}
	infile = argv[argc-2];
	outfile = argv[argc-1];

	//construct tunneled fm index
	int64_t fm_size, tfm_size;
	int64_t min_k, min_edges;
	tfm_index<> tfm;
	memory_monitor::start();
	{
		cache_config config(true, "./", util::basename(infile) );
		csa_wt<sdsl::wt_blcd<>,0xFFFFFFFF,0xFFFFFFFF> csa;
		construct( csa, infile, config, 1 );
		fm_size = size_in_bytes( csa );
		auto res = construct_tfm_index( tfm, move( csa ), config );
		tfm_size = size_in_bytes( tfm );
		min_k = res.first;
		min_edges = res.second;
	}
	memory_monitor::stop();
	
	//print infos if wanted
	if (informative) {
		//print additional information
		cout << "input_length\t" << tfm.size() << endl;
		cout << "tfm_length\t" << tfm.L.size() << endl;
		
		cout << "fm_index_size\t" << fm_size << endl;
		cout << "tfm_index_size\t" << tfm_size << endl;

		cout << "min_dbg_k\t" << min_k << endl;
		cout << "min_dbg_edges\t" << min_edges << endl;

		memory_monitor::write_memory_log<leet_format>(cout);
	}

	//serialize result
	store_to_file(tfm, outfile);
	return 0;
}

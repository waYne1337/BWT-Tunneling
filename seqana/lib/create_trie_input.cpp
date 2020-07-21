/*
 * create_trie_input.cpp for BWT Tunneling
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
#include <stdlib.h>
#include <string>
#include <tuple>

#include <sdsl/int_vector.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/io.hpp>
#include <sdsl/util.hpp>

using namespace std;
using namespace sdsl;

typedef typename sdsl::int_vector<>::size_type size_type;

void printUsage( const char *command ) {
	cerr << command << " infile outfile" << endl;
	cerr << "\tinfile: list of newline-seperated strings" << endl;
	cerr << "\toutfile: list of newline-seperated strings, where" << endl;
	cerr << "\t\t- empty lines" << endl;
	cerr << "\t\t- lines with nullbytes" << endl;
	cerr << "\t\t- strings with have another string as substring" << endl;
	cerr << "\t\t- copies of the same string except for the first occurrence" << endl;
	cerr << "\tare removed from the original list." << endl;
}

int main(int argc, char **argv ) {
	if (argc != 3) {
		printUsage(argv[0]);
		return 1;
	}
	//fetch args
	string infile = argv[1];
	string outfile = argv[2];
	csa_wt<wt_blcd<>,0x1,0xFFFFFFFF> csa;
	cache_config config(true, "./", util::basename(infile) );

	//load text
	const char* KEY_TEXT = key_text_trait<8>::KEY_TEXT;
	typedef int_vector<8> text_type;
	text_type text;
	load_vector_from_file(text, infile, 1);
	text.resize(text.size()+1);
	text[text.size()-1] = '\n'; //append an empty line at end

	//// remove empty lines and lines containing a nullbyte ///////////////
	{
		size_type linestart = 0;
		bool linehasnb = false;
		size_type j = 0;
		for (size_type i = 0; i < text.size(); i++) {
			if (text[i] == '\n') {
				if (!linehasnb && linestart != i) { //line not empty and has no nullbyte
					//copy string to front of text
					for (size_type k = linestart; k < i; k++) {
						text[j++] = text[k];
					}
					text[j++] = '\0'; //end text by a nullbyte
				}
				//prepare variables for next line
				linestart = i+1;
				linehasnb = false;
			} else if (text[i] == '\0') {
				linehasnb = true;
			}
		}
		text.resize( j );
	}

	if (text.size() == 0) {
		text.resize( 1 );
		text[0] = '\0';
	}

	//store text to cache and run csa construction
	store_to_cache(text, KEY_TEXT, config);
	register_cache_file(KEY_TEXT, config);
	construct( csa, "", config, 1 );

	//// filter strings which have another string as substring ////////////
	{
		typedef typename std::deque<std::tuple<size_type,size_type,size_type>> queue;

		//variables needed for interval_symbols function
		size_type iv_chars;
		std::vector<unsigned char> cs ( csa.sigma );
		std::vector<size_type> rank_c_i( csa.sigma );
		std::vector<size_type> rank_c_j( csa.sigma );

		//queue initialization
		queue Q;
		Q.emplace_back( 0, csa.wavelet_tree.rank( csa.size(), '\0' ), csa.size() );

		while (!Q.empty()) {
			for (size_type q = Q.size(); q > 0u; q--) {
				auto iv = Q.front();
				Q.pop_front();
			
				size_type i = get<0>(iv); //left bound
				size_type j_d = get<1>(iv); //right bound (exclusive) ending with nullbyte
				size_type j = get<2>(iv); //right bound

				//process only if entry is not marked yet
				sdsl::interval_symbols( csa.wavelet_tree, i, j_d, iv_chars, cs, rank_c_i, rank_c_j );

				if (cs.front() == '\0') { //left end of a string reached -> mark all other strings in interval
					//find leftmost string which ended (in case of copies)
					size_type lms_tpos = csa.size();
					size_type lms_sapos = csa.size();
					for (size_type rank_0 = rank_c_i.front() + 1; rank_0 <= rank_c_j.front(); rank_0++) {
						size_type s_sapos = csa.wavelet_tree.select( rank_0, '\0' );
						size_type s_tpos = csa[s_sapos];
						if (s_tpos < lms_tpos) {
							lms_tpos = s_tpos;
							lms_sapos = s_sapos;
						}
					}

					//mark all strings except of the leftmost ending one
					for (; i < j; i++) {
						if (i != lms_sapos) {
							size_type k = csa[i];
							while (text[k] != '\n' && text[k] != '\0') { //mark all positions in text to the right
								text[k++] = '\n';
							}
						}
					}
				} else { //proceed with new intervals
					for (size_type c_i = 0; c_i < iv_chars; c_i++) {
						auto c = csa.char2comp[cs[c_i]];
						size_type lb = csa.C[c] + rank_c_i[c_i];
						size_type rb_d = csa.C[c] + rank_c_j[c_i];
						size_type rb;
						if (j_d == j) {
							rb = rb_d;
						} else {
							rb = csa.C[c] + csa.wavelet_tree.rank(j, cs[c_i]);
						}
						Q.emplace_back(lb, rb_d, rb);
					}
				}					
			}
		}
	}

	//// filter strings with contained newlines ///////////////////////////
	{
		size_type linestart = 0;
		bool lineremove = false;
		size_type j = 0;
		for (size_type i = 0; i < text.size(); i++) {
			if (text[i] == '\0') {
				if (!lineremove) { //line not empty and has no nullbyte
					//copy string to front of text
					for (size_type k = linestart; k < i; k++) {
						text[j++] = text[k];
					}
					text[j++] = '\n';
				}
				//prepare variables for next line
				linestart = i+1;
				lineremove = false;
			} else if (text[i] == '\n') {
				lineremove = true;
			}
		}
		text.resize( j );
	}

	//// store result /////////////////////////////////////////////////////
	if (!store_to_plain_array<char,text_type>( text, outfile )) {
		printUsage(argv[0]);
		cerr << "Unable to write results to " << outfile << endl;
		return 1;
	}
}

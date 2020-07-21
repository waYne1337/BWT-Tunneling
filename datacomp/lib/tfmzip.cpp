/*
 * tfmzip.cpp for BWT Tunneling
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

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string.h>
#include <vector>

#include "block_compressor.hpp"
#include "twobitvector.hpp"

//from ../seqana/include/
#include "tfm_index.hpp"
#include <sdsl/io.hpp>

//post stages
#include "bcm_poststage.hpp"
#include "bw94_poststage.hpp"

using namespace std;

enum bwt_post_stage{ BW94, BCM };

//forward declarations
template<typename t_post_stage>
int tfm_compress( istream &in, ostream &out );

template<typename t_post_stage>
int tfm_decompress( istream &in, std::string &outfile );

void printUsage( char **argv ) {
	cerr << "USAGE: " << argv[0] << " [OPTIONS] INFILE OUTFILE" << endl;
	cerr << "OPTIONS:" << endl;
	cerr << "  -d\tdecompress tfm index (compression is default)." << endl;
	cerr << "    \tIf enabled, ignores all other options." << endl;
	cerr << "  -pstage [PSTAGE]\tpost stages used for the compression of the tunneled fm index. Must be one of" << endl;
	cerr << "                  \tbw94 : compression scheme from 1994 using move-to-front transform," << endl;
	cerr << "                  \t       run-length encoding and source encoding, as described" << endl;
	cerr << "                  \t       by Mike Burrows and David J. Wheeler" << endl;
	cerr << "                  \tbcm :  compression using a bwt-optimized context mixer (default)" << endl;
	cerr << "                  \t       by Ilya Muravyov" << endl;
	cerr << "INFILE:" << endl;
	cerr << "  tfm index to be compressed or decompressed if -d is set" << endl;
	cerr << "OUTFILE:" << endl;
	cerr << "  File to store the compressed (or uncompressed) data, see -d flag." << endl;
};

int main( int argc, char **argv ) {
	//analyse args
	bool compress = true; //compress or decompress
	string infile;
	string outfile;

	//check parameters
	if (argc < 3) {
		printUsage( argv );
		cerr << "At least 2 parameters expected" << endl;
		return 1;
	}

	//post bwt stage
	bwt_post_stage post_stage = BCM;

	//last option in arguments
	enum {NO, COMP, PSTAGE} last_option;
	last_option = NO;

	for (int i = 1; i < argc - 2; i++) { //analyze options
		switch (last_option) {
		case NO: //last options that require no additional parameter
		case COMP:
			if (strcmp(argv[i], "-d") == 0) { //decompress
				last_option = COMP;
				compress = false;
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
		switch (post_stage) {
		case BW94:
			fout << "w94" << endl;
			return tfm_compress<bw94_poststage>( fin, fout );
		case BCM:
			fout << "bcm" << endl;
			return tfm_compress<bcm_poststage>( fin, fout );
		}
	} else {
		fout.close();
		//read first line of input and decide what to do
		string post_stage;
		getline( fin, post_stage );
		if (post_stage == "w94") {
			return tfm_decompress<bw94_poststage>( fin, outfile );
		}
		else if (post_stage == "bcm") {
			return tfm_decompress<bcm_poststage>( fin, outfile );
		}
		else {
			cerr << "Unknown post stage " << post_stage << " used to compress index, unable to decompress" << endl;
			return 1;
		}
	}
	return 0;	
}

template<typename t_post_stage>
int tfm_compress( istream &in, ostream &out ) {
	typedef tfm_index<>::size_type size_type;
	tfm_index<> tfm;
	load( tfm, in );
	//save original text length
	block_compressor::write_primitive<size_type>( (size_type)tfm.size(), out );
	//save L
	block_compressor::write_primitive<size_type>( (size_type)tfm.L.size(), out );
	{
		std::vector<unsigned char> L( tfm.L.size() );
		for (size_type i = 0; i < L.size(); i++) {
			L[i] = (unsigned char)tfm.L[i];
		}
		t_post_stage::encode( L, out );
	}
	//save auxiliary information
	if (tfm.dout.size() != tfm.din.size()) {
		cerr << "invalid tfm index" << endl;
		return 1;
	}
	block_compressor::write_primitive<size_type>( (size_type)tfm.dout.size(), out );
	{
		twobitvector aux;
		aux.resize( tfm.dout.size() );
		for (size_type i = 0; i < aux.size(); i++) {
			aux[i] = (tfm.dout[i] << 1) + tfm.din[i];
		}
		t_post_stage::encode( aux, out );
	}
	return 0;
}

template<typename t_post_stage>
int tfm_decompress( istream &in, std::string &outfile ) {
	typedef tfm_index<>::size_type size_type;
	//prepare components
	uint64_t text_len = block_compressor::read_primitive<uint64_t>( in );
	std::string L_buf_filename = sdsl::tmp_file(outfile);
	sdsl::int_vector_buffer<8> L_buf( L_buf_filename, std::ios::out );
	sdsl::bit_vector dout;
	sdsl::bit_vector din;
	//load L
	{
		std::vector<unsigned char> L( block_compressor::read_primitive<size_type>( in ) );
		t_post_stage::decode( in, L );
		for (size_type i = 0; i < L.size(); i++) {
			L_buf[i] = L[i];
		}
	}
	//load aux
	{
		twobitvector aux;
		aux.resize( block_compressor::read_primitive<size_type>( in ) );
		t_post_stage::decode( in, aux );
		dout.resize( aux.size() );
		din.resize( aux.size() );
		for (size_type i = 0; i < aux.size(); i++) {
			dout[i] = (aux[i] >> 1);
			din[i] = (aux[i] & 1u);
		}
	}
	//construct index
	tfm_index<> tfm;
	construct_tfm_index( tfm, text_len, std::move( L_buf ), std::move( dout ), std::move( din ) );
	L_buf.close();

	//clean up and store result
	sdsl::remove( L_buf_filename );
	store_to_file(tfm, outfile);
	return 0;
}

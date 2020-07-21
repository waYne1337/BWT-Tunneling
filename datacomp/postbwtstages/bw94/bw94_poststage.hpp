/*
 * bw94_poststage.hpp for BWT Tunneling
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
 * bw94-compressor.hpp for bwt tunneling
 * Copyright (c) 2017 Uwe Baier All Rights Reserved.
 */

#ifndef BW94_POSTSTAGE_HPP
#define BW94_POSTSTAGE_HPP

#include "entropy_coder.hpp"
#include "mtf_coder.hpp"
#include "rle0_coder.hpp"

#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>

//! class which encodes a BWT with MTF + RLE0 + Entropy as second stage
class bw94_poststage {
private:
	typedef uint8_t  t_uchar_t;
	typedef uint64_t t_size_t;
	typedef uint64_t t_idx_t;
	typedef int64_t  t_bitsize_t;
public:
	//! encodes the transform t using MTF + RLE0 + Entropy
	template<class T>
	static void encode( T &t, std::ostream &out ) {
		if (t.size() == 0u)	return;
		//write alphabet
		auto alph = mtf_coder<T>::compute_alph( t );
		out.put( (t_uchar_t)alph.size() ); //store alphabet size (note that this stores 0 if full alphabet is used)
		for (t_idx_t i = 0; i < alph.size(); i++) { //and the alphabet itself
			out.put( alph[i] );
		}

		//prepare encoders
		mtf_coder<T> mtfcoder( alph );
		rle0_encoder<T> rle0coder;
		entropy_encoder<std::ostream> entcoder( out );
		entcoder.reset( alph.size() + 1 );

		for (t_idx_t i = 0; i < t.size(); ) { //do encoding
			do {
				if (i >= t.size())	break;
				//feed rle0-encoder with mtf coded input until some contents can be written
			} while (rle0coder.encode_char( mtfcoder.encode_char( t[i++] ) ));

			//move the output of the rle0coder to the entropy coder
			while (rle0coder.has_next_enc_char()) {
				entcoder.encode_char( rle0coder.next_enc_char() );
			}
		}
		entcoder.flush();
	}

	//! decodes the transform and stores it in t using MTF + RLE0 + Entropy (t must have length of output)
	template<class T>
	static void decode( std::istream &in, T &t ) {
		if (t.size() == 0u)	return;
		t_size_t alphsize = in.get();
		//check validity
		if (alphsize == 0u) {
			alphsize = std::numeric_limits<t_uchar_t>::max()+1u; //remember that on full alphabet 0 is stored
		}
		if (alphsize > t.size())
			throw std::invalid_argument("alphabet must be smaller than encoded string size");

		//read alphabet
		T alph; alph.resize( alphsize ); 
		for (t_idx_t i = 0; i < alph.size(); i++) {
			alph[i] = in.get();
		}

		//set up required decodes
		mtf_coder<T> mtfcoder( alph );
		rle0_decoder<T> rle0coder;
		entropy_decoder<std::istream> entcoder( in );
		entcoder.reset( alph.size() + 1 );

		//do decoding
		for (t_idx_t i = 0; i < t.size(); entcoder.next() ) {
			//feed rle0-decoder with input
			rle0coder.decode_char( entcoder.decode_char() );

			//fetch characters from rle0-decoder and invert mtf
			while (i < t.size() && rle0coder.has_next_char()) {
				t[i++] = mtfcoder.decode_char( rle0coder.next_char() );
			}
		}
		if (rle0coder.has_next_char()) {
			throw std::invalid_argument("encoded rle0-sequence is longer than text length");
		}
	}
};

#endif

/*
 * bwt_compressor.hpp for BWT Tunneling
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
 * bwt-compressor.hpp for bwt tunneling
 * Copyright (c) 2017 Uwe Baier All Rights Reserved.
 */

#ifndef BWT_COMPRESSOR_HPP
#define BWT_COMPRESSOR_HPP

#include "block_compressor.hpp"
#include "bwt_config.hpp"
#include "divsufsort.h"
#include "twobitvector.hpp"

#include <sdsl/util.hpp>

#include <assert.h>
#include <chrono>
#include <ios>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <stdint.h>
#include <type_traits>
#include <utility>

//! a bwt-based compressor with second stage transform as defined in t_2st_encoder
template<class tp_strategy, class t_post_stages>
class bwt_compressor : public block_compressor {
	public:
		//! constructor
		bwt_compressor() : block_compressor( t_max_size ) {};
	protected:
		virtual void compress_block( std::istream &in, std::streampos end, std::ostream &out ) const;
		virtual void decompress_block( std::istream &in, std::streampos end, std::ostream &out ) const;
};

//// COMPRESSION //////////////////////////////////////////////////////////////

template<class tp_strategy, class t_post_stages>
void bwt_compressor<tp_strategy,t_post_stages>::compress_block( std::istream &in, SDSL_UNUSED std::streampos end, std::ostream &out ) const {
	using namespace std;
	using namespace std::chrono;
	typedef typename istream::char_type schar_t;
	typedef high_resolution_clock timer;
	static_assert( is_same<
	                    typename make_unsigned<schar_t>::type,
	                    typename make_unsigned<t_uchar_t>::type
	               >::value,
	               "character types must be compatible" );

	//// GET INPUT ////////////////////////////////////////////////////////

	//get length of input
	t_size_t n = (t_size_t)(end - in.tellg());
	assert(n <= t_max_size );

	//read string from input
	t_string_t S( n );
	in.read( (schar_t *)S.data(), n );

	//// BW-TRANSFORM INPUT ///////////////////////////////////////////////

	auto start = timer::now();
	saidx_t bwt_idx = 0;
	if (bw_transform(S.data(), S.data(), NULL, (saidx_t)n, &bwt_idx) < 0) {
		throw runtime_error( string("BW Transformation failed") );
	}
	auto stop = timer::now();
	print_info("bwt_construct_time", (uint64_t)std::chrono::duration_cast<std::chrono::milliseconds>( stop - start ).count() );

	//// TUNNEL BWT ///////////////////////////////////////////////////////

	start = timer::now();
	t_idx_t tbwt_idx = bwt_idx;
	twobitvector aux;
	tp_strategy tps( S, bwt_idx );
	tps.plan();
	auto costs = tps.plan();
	print_info("num_tunnels", (uint64_t)costs.first );
	print_info("exp_tunnelcosts", (uint64_t)( costs.second / 8u) );

	auto benefit = tps.tunnel_bwt( S, aux, tbwt_idx );
	print_info("num_rle_tc", (uint64_t)benefit.first );
	print_info("exp_benefit", (uint64_t)( benefit.second / 8u) );

	stop = timer::now();
	print_info("tunneling_time", (uint64_t)std::chrono::duration_cast<std::chrono::milliseconds>( stop - start ).count() );

	//// WRITE HEADER AND ENCODING TO STREAM //////////////////////////////
	start = timer::now();

	write_primitive<t_size_t>( n, out );
	write_primitive<t_size_t>( S.size(), out );
	write_primitive<t_size_t>( aux.size(), out );
	write_primitive<t_idx_t>(  tbwt_idx, out );

	auto out_pos = out.tellp();
	t_post_stages::encode( S, out );
	print_info("size_bwt", (uint64_t)( out.tellp() - out_pos ) );
	
	out_pos = out.tellp();
	t_post_stages::encode( aux, out );
	print_info("size_aux", (uint64_t)( out.tellp() - out_pos ) );

	stop = timer::now();
	print_info("encoding_time", (uint64_t)duration_cast<milliseconds>( stop - start ).count() );
}

//// DECOMPRESSION ////////////////////////////////////////////////////////////

template<class tp_strategy, class t_post_stages>
void bwt_compressor<tp_strategy,t_post_stages>::decompress_block( std::istream &in, SDSL_UNUSED std::streampos end, std::ostream &out ) const {
	using namespace std;
	using namespace std::chrono;
	typedef typename ostream::char_type schar_t;
	typedef high_resolution_clock timer;
	static_assert( is_same<
	                    typename make_unsigned<schar_t>::type,
	                    typename make_unsigned<t_uchar_t>::type
	               >::value,
	               "character types must be compatible" );

	//// READ INPUT ///////////////////////////////////////////////////////
	auto start = timer::now();
	auto n = read_primitive<t_size_t>( in );
	auto tbwt_size = read_primitive<t_size_t>( in );
	auto aux_size = read_primitive<t_size_t>( in );
	auto tbwt_idx = read_primitive<t_idx_t>( in );
	//do some checks
	if (n > t_max_size) {
		throw invalid_argument("text(part) is too long to be decoded!");
	}
	if (tbwt_size != 0 && (tbwt_idx >= tbwt_size || tbwt_idx == 0)) {
		throw invalid_argument("invalid bwt index");
	}
	if (aux_size > tbwt_size+1) {
		throw invalid_argument("aux size is longer than tbwt size");
	}
	t_string_t tbwt; tbwt.resize( tbwt_size );
	twobitvector aux; aux.resize( aux_size );
	t_post_stages::decode( in, tbwt );
	t_post_stages::decode( in, aux );
	auto stop = timer::now();
	print_info("decoding_time", (uint64_t)duration_cast<milliseconds>( stop - start ).count() );

	//// INVERT TUNNELED BWT //////////////////////////////////////////////
	start = timer::now();
	tp_strategy::invert_tbwt( std::move(tbwt), std::move(aux), n, tbwt_idx, out );
	stop = timer::now();
	print_info("inversion_time", (uint64_t)duration_cast<milliseconds>( stop - start ).count() );
}

#endif

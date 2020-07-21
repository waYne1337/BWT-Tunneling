/*
 * tp_strategy_none.hpp for BWT Tunneling
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

#ifndef TP_STRATEGY_NONE_HPP
#define TP_STRATEGY_NONE_HPP

#include "bwt_config.hpp"
#include "divsufsort.h"
#include "twobitvector.hpp"

#include <utility>

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

class tp_strategy_none {
public:
	tp_strategy_none( SDSL_UNUSED const t_string_t &L, SDSL_UNUSED t_idx_t bwt_idx ) {};

	std::pair<t_size_t,t_bitsize_t> plan() {
		return std::pair<t_size_t,t_bitsize_t>( 0u, 0u );
	};

	std::pair<t_size_t,t_bitsize_t> tunnel_bwt( SDSL_UNUSED t_string_t &bwt, SDSL_UNUSED twobitvector &aux, SDSL_UNUSED t_idx_t &bwt_idx ) {
		return std::pair<t_size_t,t_bitsize_t>( 0u, 0u );
	};

	static void invert_tbwt( t_string_t &&tbwt, twobitvector &&aux, t_size_t n,
                                 t_idx_t tbwt_idx, std::ostream &out );
};

//// INVERTING A TUNNELED BWT /////////////////////////////////////////////////

void tp_strategy_none::invert_tbwt( t_string_t &&tbwt, SDSL_UNUSED twobitvector &&aux, SDSL_UNUSED t_size_t n,
                                           t_idx_t tbwt_idx, std::ostream &out ) {
	typedef typename std::ostream::char_type schar_t;
	static_assert( std::is_same<
	                    typename std::make_unsigned<schar_t>::type,
	                    typename std::make_unsigned<t_uchar_t>::type
	               >::value,
	               "character types must be compatible" );

	if (tbwt.size() != 0 && (tbwt_idx >= tbwt.size() || tbwt_idx == 0)) {
		throw std::invalid_argument("tbwt index is invalid");
	}
	if (inverse_bw_transform(tbwt.data(), tbwt.data(), NULL,
	                         (saidx_t)tbwt.size(), (saidx_t)tbwt_idx) < 0) {
		throw std::invalid_argument( "Inverse BW Transformation failed" );		
	}
	out.write( (const schar_t *)tbwt.data(), tbwt.size() );
}

#endif

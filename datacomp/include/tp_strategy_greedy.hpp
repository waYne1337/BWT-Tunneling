/*
 * tp_strategy_greedy.hpp for BWT Tunneling
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

#ifndef TP_STRATEGY_GREEDY_HPP
#define TP_STRATEGY_GREEDY_HPP

#include "bwt_config.hpp"
#include "tp_strategy_lmrtpi.hpp"

#include <sdsl/bits.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

class tp_strategy_greedy : public tp_strategy_lmrtpi {
public:
		tp_strategy_greedy( const t_string_t &L, t_idx_t bwt_idx ) : tp_strategy_lmrtpi( L, bwt_idx ) {
	};

	virtual std::pair<t_size_t,t_bitsize_t> plan() {
		//compute rating
		std::vector<t_size_t> RPTC;
		compute_rating( RPTC );

		//set up a list of run identifiers
		sdsl::int_vector<> SPI( r, 0, sdsl::bits::hi( r ) + 1 );
		sdsl::util::set_to_id( SPI );

		//sort the identifiers according to the rating array (descending s.t. best ratings are in front)
		sort( SPI.begin(), SPI.end(), [&]( t_idx_t i, t_idx_t j ) { return RPTC[i] > RPTC[j]; } );

		//find the best value for t
		t_size_t t_opt = 0;
		t_bitsize_t ben_opt = benefit(0u) - benefit(0u);
		for (t_size_t t = 1, tc = 0; t <= r; t++) {
			t_idx_t k = SPI[t-1];
			tc += RPTC[k];
			t_bitsize_t ben = benefit(tc) - cost(t);
			if (ben > ben_opt) {
				t_opt = t;
				ben_opt = ben;
			}
		}

		//clear prefix intervals not belonging to best choice
		for (t_size_t t = t_opt; t < r; t++) {
			t_idx_t k = SPI[t];
			RPE[k] = run_lf.lfr(k);
		}
		return std::pair<t_size_t,t_bitsize_t>( t_opt, cost(t_opt) );
	};			
};

#endif

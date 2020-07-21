/*
 * tp_strategy_hirsch.hpp for BWT Tunneling
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

#ifndef TP_STRATEGY_HIRSCH_HPP
#define TP_STRATEGY_HIRSCH_HPP

#include "bwt_config.hpp"
#include "tp_strategy_lmrtpi.hpp"

#include <sdsl/bits.hpp>

#include <algorithm>
#include <math.h>
#include <vector>

class tp_strategy_hirsch : public tp_strategy_lmrtpi {
private:
	t_size_t compute_MT( t_size_t tc ) {
		t_bitsize_t p = std::max( 0.0, (double)( ( tc * log2_2nrle_rc - 2) / 4 ));
		p = std::min( (t_bitsize_t)sdsl::bits::hi( rhg1 ) + 1, p ); //limit amount of left-shifting
		return (t_size_t)(((rhg1 + 1) / ((1u << p) + 2)) - 0.5);
	};
public:
	tp_strategy_hirsch( const t_string_t &L, t_idx_t bwt_idx ) : tp_strategy_lmrtpi( L, bwt_idx ) {
	};

	virtual std::pair<t_size_t,t_bitsize_t> plan() {
		//compute rating
		std::vector<t_size_t> RPTC;
		compute_rating( RPTC );

		//replace ratings with the amount of tunnels needed such that it is worth
		//to tunnel the prefix interval
		for (t_idx_t k = 0; k < r; k++) {
			RPTC[k] = compute_MT( RPTC[k] );
		}
		//set up a counting array and count values in RPTC
		std::vector<t_size_t> C( compute_MT( 0u ), 0 );
		for (t_idx_t k = 0; k < r; k++) {
			if (RPTC[k] < C.size()) {
				++C[RPTC[k]];
			}
		}
		//find largest value t_opt such that t_opt prefix intervals have more benefit than cost
		t_size_t t_opt = 0;
		for (t_size_t t = 1; t < C.size(); t++) {
			C[t] += C[t-1];
			if (C[t] >= t) {
				t_opt = t;
			}
		}
		//clear prefix intervals which have less benefit than cost
		for (t_idx_t k = 0; k < r; k++) {
			if (RPTC[k] > t_opt) {
				RPE[k] = run_lf.lfr(k);
			}
		}
		return std::pair<t_size_t,t_bitsize_t>( C[t_opt], cost(t_opt) );
	};
};

#endif

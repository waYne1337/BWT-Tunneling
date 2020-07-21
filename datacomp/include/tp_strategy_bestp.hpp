/*
 * tp_strategy_bestp.hpp for BWT Tunneling
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

#ifndef TP_STRATEGY_BESTP_HPP
#define TP_STRATEGY_BESTP_HPP

#include "bwt_config.hpp"
#include "tp_strategy_lmrtpi.hpp"

class tp_strategy_bestp : public tp_strategy_lmrtpi {
public:
	static t_size_t percentage; //percentage to be used

	tp_strategy_bestp( const t_string_t &L, t_idx_t bwt_idx ) : tp_strategy_lmrtpi( L, bwt_idx ) {
	};

	virtual std::pair<t_size_t,t_bitsize_t> plan() {
		//compute rating
		std::vector<t_size_t> RPTC;
		compute_rating( RPTC );

		//set up a list of run identifiers
		sdsl::int_vector<> SPI( r, 0, sdsl::bits::hi( r ) + 1 );
		sdsl::util::set_to_id( SPI );

		//count amount of tunnelable prefix intervals
		t_size_t num_tunnels = 0;
		for (t_size_t k = 0; k < r; k++) {
			if (RPTC[k] != 0)	++num_tunnels;
		}

		//sort the identifiers according to the rating array (descending s.t. best ratings are in front)
		sort( SPI.begin(), SPI.end(), [&]( t_idx_t i, t_idx_t j ) { return RPTC[i] > RPTC[j]; } );

		//clear prefix intervals not belonging to the best percentage
		num_tunnels = ((t_bitsize_t)num_tunnels * (t_bitsize_t)percentage) / 100u;
		for (t_size_t t = num_tunnels; t < r; t++) {
			t_idx_t k = SPI[t];
			RPE[k] = run_lf.lfr(k);
		}
		return std::pair<t_size_t,t_bitsize_t>( num_tunnels, cost(num_tunnels) );
	};
};

t_size_t tp_strategy_bestp::percentage = 10;

#endif

/*
 * tp_strategy_greedy_update.hpp for BWT Tunneling
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

#ifndef TP_STRATEGY_GREEDY_UPDATE_HPP
#define TP_STRATEGY_GREEDY_UPDATE_HPP

#include <sdsl/bits.hpp>
#include <sdsl/int_vector.hpp>

#include "bwt_config.hpp"
#include "lheap.hpp"
#include "tp_strategy_lmrtpi.hpp"
#include "twobitvector.hpp"

class tp_strategy_greedy_update : public tp_strategy_lmrtpi {
private:
	//! utility function to execute some command for all runs
	//! surrounding the columns of a prefix interval except of the first one
	template<class column_handler>
	void enumerate_columns( t_idx_t k, column_handler &handler ) {
		
		auto x = run_lf.lfr(k);
		while (x != RPE[k]) {
			auto k_ = run_lf.run_of(x);
			handler( k, k_ );
			x = run_lf.lfr( k_ ) + (x - run_lf.start( k_ ) );
		}
	};
		
public:
	tp_strategy_greedy_update( const t_string_t &L, t_idx_t bwt_idx ) : tp_strategy_lmrtpi( L, bwt_idx ) {
	};

	virtual std::pair<t_size_t,t_bitsize_t> plan() {
		//compute rating
		std::vector<t_size_t> RPTC;
		compute_rating( RPTC );

		//set up a map of prefix intervals pointing through runs
		sdsl::int_vector<> ov_graph( run_lf.n, 0, sdsl::bits::hi( r ) + 1 );
		//and an array storing the length of the prefix interval
		sdsl::int_vector<> length( r, 0, sdsl::bits::hi( run_lf.n ) + 1 );

		//initialize both arrays
		auto col_handler = [&]( t_idx_t k, t_idx_t k_ ) {
			length[k]++;
			if (run_lf.height(k_) > run_lf.height(k)) {
				//add k t the overlap entries of k_
				auto i = run_lf.start( k_ );
				auto off = ++ov_graph[i]; //start position stores number of entries
				ov_graph[i + off] = k; //entries are stored below
			}
		};
		for (t_idx_t k = 0; k < r; k++) {
			enumerate_columns( k, col_handler );
		}

		//initialize a heap containing all elements
		std::vector<t_idx_t> H; H.reserve( r );
		twobitvector HS;	HS.resize( r );
		for (t_idx_t k = 0; k < r; k++) {
			if (RPTC[k] > 0u && length[k] > 2u) {
				H.push_back( k );
				HS[k] = lheap_vstate::unchanged;
			} else {
				HS[k] = lheap_vstate::empty;
			}
		}
		auto pi_cmp = [&]( t_idx_t k1, t_idx_t k2 ) { //rating comparison
			return RPTC[k1] < RPTC[k2];
		};
		auto pi_state = [&]( t_idx_t k ) { //state indication in heap
			return HS[k];
		};
		make_lheap( H.begin(), H.end(), pi_cmp );

		//greedily choose biggest prefix intervals and update overlapping prefix intervals
		t_size_t t = 0;
		t_size_t tc = 0;
		t_size_t t_opt = 0;
		t_bitsize_t ben_opt = 0;
		auto SPI = H.rbegin(); //array to store sorted prefix intervals
		for (auto e = H.end(); e != H.begin(); ) {
			t_idx_t k = H.front(); //get prefix interval with maximal score

			//check for a new best choice
			tc += RPTC[k];
			++t;
			t_bitsize_t ben = benefit(tc) - cost(t);
			if (ben > ben_opt) {
				t_opt = t;
				ben_opt = ben;
			}

			//update ratings of overlapping prefix intervals where k points through
			auto update_inner = [&]( t_idx_t k_cur, t_idx_t k_inner ) {
				if (HS[k_inner] != lheap_vstate::empty) {
					t_size_t rptc_diff = ((t_bitsize_t)(length[k_inner]-2)) * RPTC[k_cur] / (length[k_cur] - 2);
					if (rptc_diff < RPTC[k_inner]) {
						RPTC[k_inner] -= rptc_diff;
						HS[k_inner] = lheap_vstate::decreased;
					} else {
						RPTC[k_inner] = 0;
						HS[k_inner] = lheap_vstate::empty;
					}
				}
			};
			enumerate_columns( k, update_inner );

			//update ratings of overlapping prefix intervals pointing through k
			auto update_outer = [&]( t_idx_t k_cur, t_idx_t k_outer ) {
				if (HS[k_outer] != lheap_vstate::empty) {
					t_bitsize_t new_length = length[k_outer] - ( length[k_cur] - 2 );
					t_size_t new_rptc = ( (new_length-2) * RPTC[k_outer] ) / ( length[k_outer] - 2 );
					if (new_length > 2u && new_rptc > 0u) {
						length[k_outer] = new_length;
						RPTC[k_outer] = new_rptc;
						HS[k_outer] = lheap_vstate::decreased;
					} else {
						RPTC[k_outer] = 0;
						HS[k_outer] = lheap_vstate::empty;
					}
				}
			};
			//list of prefix intervals can be found in ov_graph[run_lf.start(k)+1,run_lf.start(k)+ov_graph[run_lf.start(k)]]
			t_idx_t j = run_lf.start(k) + ov_graph[run_lf.start(k)];
			for (t_idx_t i = run_lf.start(k) + 1; i <= j; i++) {
				t_idx_t k_outer = ov_graph[i];
				update_outer( k, k_outer );
			}

			//remove prefix interval from heap, and store it in SPI (similar to heapsort)
			e = pop_lheap_nomove(H.begin(), e, pi_state, pi_cmp);
			HS[k] = lheap_vstate::empty;
			*(SPI++) = k;
		}

		//set the heap state of best choice to unchanged
		SPI = H.rbegin();
		for (t_idx_t t = 0; t < t_opt; t++) {
			auto k = SPI[t];
			HS[k] = lheap_vstate::unchanged;
		}

		//remove markings of all prefix intervals which should not be tunneled
		for (t_idx_t k = 0; k < r; k++) {
			if (HS[k] != lheap_vstate::unchanged) {
				RPE[k] = run_lf.lfr( k );
			}
		}
		return std::pair<t_size_t,t_bitsize_t>( t_opt, cost(t_opt) );
	};
};

#endif

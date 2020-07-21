/*
 * run_lf_support.hpp for BWT Tunneling
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

#ifndef RUN_LF_SUPPORT_HPP
#define RUN_LF_SUPPORT_HPP

#include <algorithm>
#include <limits>
#include <vector>

#include "bwt_config.hpp"

//! support structure for bwt navigation and bwt run support
class run_lf_support {
	private:
		t_size_t m_runs; //number of logical runs
		t_size_t m_idx_runs; //number of runs (indexed BWT)
		t_idx_t m_bwt_idx; //bwt index
		t_size_t m_n; //logical text length
		t_size_t m_idx_n; //text length (indexed BWT)
		t_size_t m_sigma; //size of alphabet
		t_size_t m_max_char_val; //maximal value of an element in alphabet

		std::vector<t_idx_t> m_lfr; //lf, only for the start of runs
		std::vector<t_idx_t> m_rs; //start positions of all runs, sorted ascending.
		                           //additionally, m_rs[m_runs] = n+1 holds.

	public:
		//! constructor, expects a indexed BWT and its primary index.
		run_lf_support( const t_uchar_t *bwt, t_size_t _n, t_idx_t _idx );

		//! logical number of runs in BWT
		const t_size_t &runs = m_runs;

		//! number of runs in the indexed BWT
		const t_size_t &idx_runs = m_idx_runs;

		//! primary index of the bwt
		const t_idx_t &bwt_idx = m_bwt_idx;

		//! logical length of text
		const t_size_t &n = m_n;

		//! real text length (also length of indexed BWT)
		const t_size_t &idx_n = m_idx_n;

		//! size of alphabet in text
		const t_size_t &sigma = m_sigma;

		//! maximal value of an element in alphabet
		const t_size_t &max_char_val = m_max_char_val;

		//! utility function, returns lf at the start of the given run
		t_idx_t lfr( t_idx_t r ) const {
			return m_lfr[r];
		};

		//! utility function, returns the start of a run
		t_idx_t start( t_idx_t r ) const {
			return m_rs[r];
		};

		//! function returns the run to which position i belongs,
		//  or a value >= runs if i does not belong to any run (e.g. i < 0 or i >= n)
		t_idx_t run_of( t_idx_t i ) const;

		//! utility function, computes height of a run
		t_size_t height( t_idx_t r ) const {
			return m_rs[r+1]-m_rs[r];
		};

		//! utility function, computes exclusive end of a run
		t_idx_t end( t_idx_t r ) const {
			return m_rs[r+1];
		};

		//! utility function, converts a position in the indexed bwt to
		//! a position in the logical bwt
		t_idx_t idx_to_log( t_idx_t p_idx ) const {
			return (p_idx < bwt_idx) ? p_idx : p_idx + 1;
		};

		//! utility function, converts a logical position in the bwt to
		//! a position in the indexed bwt
		t_idx_t log_to_idx( t_idx_t p_log ) const {
			return (p_log <= bwt_idx) ? p_log : p_log - 1;
		};
};

t_idx_t run_lf_support::run_of( t_idx_t i ) const {
	//use binary search with runstart - array
	auto it = std::upper_bound( m_rs.begin(), m_rs.end(), i );
	return (t_idx_t)(it - m_rs.begin()) - 1;
}

run_lf_support::run_lf_support( const t_uchar_t *bwt, t_size_t _n, t_idx_t idx ) {
	//init some basic variables
	m_bwt_idx = idx;
	m_idx_n = _n;
	m_idx_runs = 0;
	m_sigma = 0;
	m_max_char_val = 0;

	//build C Array and count runs
	std::vector<t_size_t> C( std::numeric_limits<t_uchar_t>::max() + 1 );
	t_idx_t borders[] = {bwt_idx,idx_n};
	t_idx_t i = 0;
	for (t_idx_t b : borders) { //to split runs at primary index
		t_uchar_t lastchar = (i < idx_n) ? bwt[i]+1 : 0;
		while (i < b) {
			if (lastchar != bwt[i]) { //start of a run
				++m_idx_runs;
				lastchar = bwt[i];
			}
			++C[lastchar];
			++i;
		}
	}
	m_n = idx_n + 1;
	m_runs = idx_runs + 1; //for bwt index

	//build cumulative sums of the C array
	t_idx_t l = 1; //for bwt index
	for (t_idx_t c = 0; c < C.size(); c++) {
		auto tmp = C[c];
		C[c] = l;
		l += tmp;
		if (tmp > 0) {
			++m_sigma;
			m_max_char_val = c;
		}
	}

	//compute LF
	m_lfr.reserve( m_runs + 1 );
	m_rs.reserve( m_runs + 1 );
	i = 0;
	t_idx_t i_log = 0; //logical position of i
	for (t_idx_t b : borders) { //to split runs at primary index
		t_uchar_t lastchar = (i < n) ? bwt[i]+1 : 0;
		while (i < b) {
			if (lastchar != bwt[i]) { //start of a run
				m_rs.push_back( i_log ); //store start of run

				lastchar = bwt[i];
				m_lfr.push_back( C[lastchar] );
			}
			++C[lastchar];
			++i; ++i_log;
		}
		//add a terminator to both lfr and rs (for both primary index and n)
		m_rs.push_back( i_log++ );
		m_lfr.push_back( 0 );
	}
}

#endif

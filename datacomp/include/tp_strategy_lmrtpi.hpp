/*
 * tp_strategy_rtpi.hpp for BWT Tunneling
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

#ifndef TP_STRATEGY_LMRTPI_HPP
#define TP_STRATEGY_LMRTPI_HPP

#include <limits>
#include <math.h>
#include <ostream>
#include <stack>
#include <stdexcept>
#include <utility>
#include <vector>

#include <sdsl/bits.hpp>

#include "aux_encoding.hpp"
#include "bwt_config.hpp"
#include "run_lf_support.hpp"
#include "twobitvector.hpp"

//! tunneling strategy considering length-maximal run-terminated prefix intervals
class tp_strategy_lmrtpi {
private:
	void compute_lmrtpis();

	static void transform_aux( const t_string_t &tbwt, twobitvector &aux, t_idx_t tbwt_idx );
	static void retransform_aux( const t_string_t &tbwt, twobitvector &aux, t_idx_t tbwt_idx );

	static t_size_t rle_len( const t_string_t &L );
protected:
	//run-lf support
	run_lf_support run_lf;
	
	//variables
	t_size_t n_rle;//length of run-length encoding
	t_size_t r;    //number of runs
	t_size_t rhg1; //number of runs with height  > 1
	t_size_t rc;   //overall amount of run characters in BWT
	double log2_2nrle_rc; //log_2( (2*n_rle) / rc )

	//prefix interval array
	std::vector<t_size_t> RPE;

	//computes the rating array
	void compute_rating(std::vector<t_size_t> &RPTC);

public:
	//! constructor, get information
	tp_strategy_lmrtpi( const t_string_t &L, t_idx_t bwt_idx ) : run_lf( (const t_uchar_t *)L.data(), L.size(), bwt_idx ) {
		//compute variables
		r = run_lf.runs;
		rhg1 = 0;
		rc = 0;
		n_rle = rle_len( L );
		rc = n_rle - r;
		for (t_idx_t i = 0; i < r; i++) {
			if (run_lf.height(i) > 1u) {
				++rhg1;
			}
		}
		log2_2nrle_rc = 1.0 + log1p( r / (double)rc ) / log( 2 );

		//create RPE array
		compute_lmrtpis();
	};

	//! benefit function for tc removed characters
	t_bitsize_t benefit( t_size_t tc ) const {
		return round( tc * log2_2nrle_rc );
	};

	//! cost function for t tunnels
	t_bitsize_t cost( t_size_t t ) const {
		return (rhg1 + 1 >= 4*t + 2)
		       ? (t + 0.5) * (6 + 4 * log2( (rhg1 + 1) / (double)( 2*t + 1) - 1.0 ))
		       : (t + 0.5) * 6;
	};

	//! planning, returns number of prefix intervals to be tunneled and expected cost
	virtual std::pair<t_size_t,t_bitsize_t> plan() = 0;

	//! tunneling according to the plan. Returns the number of removed characters from the RLE representation
	//! and the expected benefit in bits
	std::pair<t_size_t,t_bitsize_t> tunnel_bwt( t_string_t &bwt, twobitvector &aux, t_idx_t &tbwt_idx );

	//! invert a tunneled BWT
	static void invert_tbwt( t_string_t &&tbwt, twobitvector &&aux, t_size_t n,
                                 t_idx_t tbwt_idx, std::ostream &out );
};

//// COMPUTATION OF THE LENGTH OF A RUN-LENGTH ENCODING ///////////////////////
t_size_t tp_strategy_lmrtpi::rle_len( const t_string_t &L ) {
	if (L.size() == 0)	return 0;

	t_uchar_t rchar = L[0] + 1; //character of current run
	t_size_t rlen = 1; //length of current run
	t_size_t rle_len = 0; //overall length
	for (t_idx_t i = 0; i < L.size(); i++) {
		if (rchar != L[i]) { //end of a run
			rle_len += 1 + sdsl::bits::hi( rlen );
			rchar = L[i];
			rlen = 1;
		} else {
			++rlen;
		}
	}
	rle_len += sdsl::bits::hi( rlen ); //add characters of last run (1 is added for the zero-th run)
	return rle_len;
}

//// COMPUTATION OF LENGTH-MAXIMAL RUN-TERMINATED PREFIX INTERVALS ////////////
void tp_strategy_lmrtpi::compute_lmrtpis() {
	//create helpful arrays
	RPE.resize( r );
	std::vector<t_idx_t> PE( r ); //prefix interval end

	//initialize arrays and a stack
	for (t_idx_t k = 0; k < r; k++) {
		PE[k] = run_lf.lfr(k);
		RPE[k] = run_lf.lfr(k);
	}
	std::stack<t_idx_t> s;

	//run algorithm
	for (t_idx_t k = 0; k < r; k++) {
		if (run_lf.height(k) >= 2u) {
			s.push( k );
			do {
				auto r_ = run_lf.run_of( PE[s.top()] );
				if (PE[s.top()] + run_lf.height(s.top()) <= run_lf.end(r_)) {
					//stacktop prefix interval can be extended
					s.push( r_ );
				} else {
					//prefix interval on stacktop is maximal
					r_ = s.top();
					s.pop();

					if (!s.empty()) {
						//adapt end of new stacktop prefix interval
						PE[s.top()] = PE[r_] + (PE[s.top()] - run_lf.start(r_));

						//adapt end of prefix intervals
						if (run_lf.height(r_) == run_lf.height(s.top())) {
							RPE[s.top()] = RPE[r_];
							RPE[r_] = run_lf.lfr(r_);
						}
					}
				}
			} while (!s.empty());
		}
	}
};

//// COMPUTATION OF RATING ARRAY //////////////////////////////////////////////
void tp_strategy_lmrtpi::compute_rating(std::vector<t_size_t> &RPTC) {
	RPTC.resize( r );
	for (t_idx_t k = 0; k < r; k++) {
		RPTC[k] = 0;
		if (RPE[k] != run_lf.lfr(k)) {	
			auto h = run_lf.height(k);
			auto x = run_lf.lfr(k);
			while (x != RPE[k]) {
				auto r_ = run_lf.run_of(x);
				auto h_ = run_lf.height(r_);
				RPTC[k] += sdsl::bits::hi( h_ ) - sdsl::bits::hi( h_ - h + 1 );
				x = run_lf.lfr( r_ ) + (x - run_lf.start( r_ ) );
			}
			RPTC[k] -= sdsl::bits::hi( h );
		}
	}
};

//// TRANSFORM AUX ////////////////////////////////////////////////////////////
void tp_strategy_lmrtpi::transform_aux( const t_string_t &tbwt, twobitvector &aux, t_idx_t tbwt_idx ) {
	if (tbwt.size() == 0)	return;

	//transfer aux to a run-based representation
	t_idx_t j = 0;
	std::array<t_idx_t,3> bounds = {(t_idx_t)0u,
		                   tbwt_idx,
		                   (t_idx_t)tbwt.size()}; //don't forget primary index run
	for (t_idx_t ib = 0; ib != 2; ib++) {
		t_idx_t i = bounds[ib];
		t_idx_t e = bounds[ib+1];

		bool newrun = true;
		while (++i < e) {
			if (tbwt[i] != tbwt[i-1]) { //new run detected
				newrun = true;
			}
			else if (newrun) { //first run-character of new run
				aux[j++] = aux[i]; //copy aux-value of runs with height > 1
				newrun = false;
			}
		}
	}
	aux.resize( j );
}

//// RETRANSFORM AUX //////////////////////////////////////////////////////////
void tp_strategy_lmrtpi::retransform_aux( const t_string_t &tbwt, twobitvector &aux, t_idx_t tbwt_idx ) {
	if (tbwt.size() == 0)	return;

	//decode run-based aux
	std::array<t_idx_t,3> bounds = {
		(t_idx_t)tbwt.size(),
		tbwt_idx,
		(t_idx_t)0u }; //don't forget primary index run

	t_idx_t j = aux.size();
	aux.resize( tbwt.size() + 1 );
	aux[tbwt.size()] = aux_encoding::REG;

	for (t_idx_t ib = 0; ib != 2; ib++) {
		t_idx_t i = bounds[ib];
		t_idx_t e = bounds[ib+1];

		bool newrun = true;
		while (--i > e) {
			if (tbwt[i] != tbwt[i-1]) { //start of a run
				aux[i] = aux_encoding::REG;
				newrun = true;
			}
			else {
				if (newrun) {
					if (j-- == 0u) throw std::invalid_argument("invalid aux encoding");
					newrun = false;
				}
				aux[i] = aux[j];
			}
		}
		aux[i] = aux_encoding::REG; //set flag for start of last run
	}
}

//// TUNNEL A BWT /////////////////////////////////////////////////////////////
std::pair<t_size_t,t_bitsize_t> tp_strategy_lmrtpi::tunnel_bwt( t_string_t &bwt, twobitvector &aux, t_idx_t &tbwt_idx ) {

	//resize auxiliary bit vector to cover enough space
	aux.resize( run_lf.idx_n+1 );
	for (t_idx_t i = 0; i < aux.size(); i++) {
		aux[i] = aux_encoding::REG;
	}			

	//mark each tunnel in auxiliary structure
	std::vector<t_idx_t> intervals;
	for (t_idx_t k = 0; k < r; k++) {
		if (RPE[k] != run_lf.lfr(k)) { //tunnel prefix interval
			//save which rows of prefix interval were not tunneled yet
			intervals.clear();
			auto lastaux = aux_encoding::REM;
			for (t_idx_t i = run_lf.log_to_idx(run_lf.start(k)+1);
				     i < run_lf.log_to_idx(run_lf.end(k)); i++) {
				if (aux[i] != lastaux) {
					lastaux = aux[i];
					intervals.push_back( i - run_lf.log_to_idx(run_lf.start(k)) );
				}
			}
			intervals.push_back( run_lf.height(k) );

			//prepare marking in aux
			auto cur = run_lf.lfr(k);    //current run
			auto last = run_lf.start(k); //previous run (seen in text order)
			while (cur != RPE[k]) {
				//clear cntL for all intervals in previous run
				for (t_idx_t i = 1; i < intervals.size(); i += 2 ) {
					t_idx_t i_s = run_lf.log_to_idx(last + intervals[i-1]); //interval start
					t_idx_t i_e = run_lf.log_to_idx(last + intervals[i]  ); //interval end
					do {
						aux[i_s] = aux[i_s] | aux_encoding::IGN_L;
					} while (++i_s < i_e);
				}

				//move on cur and last by 1 column
				last = cur;
				t_idx_t cur_r = run_lf.run_of( cur );
				cur = (aux[run_lf.log_to_idx( cur + intervals.front() )] == aux_encoding::IGN_L) //prefix interval of run cur_r was tunneled already
				    ? RPE[cur_r]        + (cur - run_lf.start(cur_r))  //jump over prefix interval if tunneled already
				    : run_lf.lfr(cur_r) + (cur - run_lf.start(cur_r)); //otherwise proceed stepwise

				//clear cntF for all intervals in current run (before both pointers were moved)
				for (t_idx_t i = 1; i < intervals.size(); i += 2 ) {
					t_idx_t i_s = run_lf.log_to_idx(last + intervals[i-1]);
					t_idx_t i_e = run_lf.log_to_idx(last + intervals[i]  );
					do {
						aux[i_s] = aux[i_s] | aux_encoding::SKP_F;
					} while (++i_s < i_e);
				}
			}
			//set end of prefix interval k one position to right (i.e. one application of inverse LF),
			// such that prefix interval jumping as shown above works correct
			RPE[k] = last;
		}
	}

	//remove doubly marked entries from BWT and auxiliary structure
	t_idx_t borders[] = { run_lf.bwt_idx, run_lf.idx_n };
	t_idx_t p = 0; //position in bwt and auxiliary data structure
	t_idx_t i = 0; //position in original bwt
	for (auto b : borders) {
		tbwt_idx = p; //also, compute new position of primary index
		while (i < b) {
			if (aux[i] != aux_encoding::REM) { //copy entries which won't be removed
				bwt[p] = bwt[i];
				aux[p++] = aux[i];
			}
			++i;
		}
	}
	//trim both bwt and aux to correct sizes and add a terminator to aux
	bwt.resize( p );
	aux[p++] = aux_encoding::REG;	aux.resize( p );

	transform_aux( bwt, aux, tbwt_idx );

	//measure removed characters from RLE encoding
	t_size_t tc = n_rle - rle_len( bwt );
	return std::pair<t_size_t,t_bitsize_t>( tc, benefit(tc) );
}

//// INVERTING A TUNNELED BWT /////////////////////////////////////////////////

void tp_strategy_lmrtpi::invert_tbwt( t_string_t &&tbwt, twobitvector &&aux, t_size_t n,
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

	retransform_aux( tbwt, aux, tbwt_idx );

	//count character frequencies
	std::vector<t_size_t> C( std::numeric_limits<t_uchar_t>::max() + 1 );
	for (t_idx_t i = 0; i < tbwt.size(); i++) {
		if (aux[i] != aux_encoding::IGN_L) { //skip entries which are ignored
			++C[tbwt[i]];
		}
	}
	//compute start positions
	t_size_t j = 0;
	for (t_size_t i = 0; i < C.size(); i++) {
		auto cnt = C[i];
		C[i] = j;
		while (cnt > 0) {
			if (++j >= aux.size()) {
				throw std::invalid_argument("auxiliary structure is invalid");
			}				
			if (aux[j] != aux_encoding::SKP_F) { //skip empty positions
				--cnt;
			}
		}
	}

	//NOTE: the following code requires that aux[tbwt.size()] == aux_encoding::REG
	if (aux[tbwt.size()] != aux_encoding::REG) {
		throw std::invalid_argument("auxiliary structure is invalid");
	}

	//// INVERTITION USING PHI ////////////////////////////////////////////
	
	//compute PHI
	std::vector<t_idx_t> PHI( tbwt.size() );
	for (t_idx_t i = 0; i < tbwt.size(); i++) {
		if (aux[i] != aux_encoding::IGN_L) {
			j = C[tbwt[i]];
			if (j < tbwt_idx) {
				//skip empty positions
				for (t_idx_t k = 1; aux[++j] == aux_encoding::SKP_F; k++) {
					PHI[j] = k; //save distance to previous regular entry
				}
				//set PHI
				if (j >= tbwt.size()) {
					throw std::invalid_argument("auxiliary structure is invalid");
				}
				if (j < tbwt_idx)	PHI[j] = i;
				else             	PHI[0] = i; //save start index
			} else {
				//set PHI
				if (j >= tbwt.size()) {
					throw std::invalid_argument("auxiliary structure is invalid");
				}
				PHI[j] = i;
				//skip empty positions
				for (t_idx_t k = 1; aux[++j] == aux_encoding::SKP_F; k++) {
					PHI[j] = k; //save distance to previous regular entry
				}
			}
			C[tbwt[i]] = j;
		}
	}
	//invert tunneled bwt using a stack
	std::stack<t_idx_t> stck;
	j = 0; //start at saved start index
	for (t_idx_t i = 0; i < n; i++) {
		j = PHI[j];
		out.put( (schar_t)tbwt[j] );
		if ( aux[j+1] == aux_encoding::IGN_L ) { //end of a tunnel (reverse order)
			if (stck.empty()) {
				throw std::invalid_argument("missing end of a tunnel");
			}
			else {
				j += stck.top();
				stck.pop();
			}
		}
		else if ( aux[j] == aux_encoding::SKP_F ) { //start of a tunnel (reverse order)
			stck.push( PHI[j] ); //save distance to uppermost row of block
			j -= PHI[j];
		}
		else if ( aux[j+1] == aux_encoding::SKP_F ) { //start of a tunnel, being at the uppermost row (reverse order)
			stck.push( 0 );
		}
	}
	if (!stck.empty()) {
		throw std::invalid_argument("missing start of a tunnel");
	}
}

#endif

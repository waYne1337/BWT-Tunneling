/*
 * trie_construct_algorithms.hpp for BWT Tunneling
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
//ALGORITHMS COME FROM:
/*
 * xbwt_construct_algorithm.hpp for trickier xbwt tricks
 * Copyright (c) 2018 Stefan Stauß, Uwe Baier All Rights Reserved.
 *
 * BWT Tunneling
 * Copyright (c) 2020 Uwe Baier All Rights Reserved.
 */

#ifndef TRIE_CONSTRUCT_ALGORITHMS_HPP
#define TRIE_CONSTRUCT_ALGORITHMS_HPP

#include <algorithm>
#include <utility>
#include <type_traits>
#include <queue>
#include <utility>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/construct_config.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/util.hpp>

#include "dbg_algorithms.hpp"
#include "succinct_counter.hpp"
#include "tries.hpp"
#include "trie_xbwt.hpp"
#include "trie_txbwt.hpp"

//forward declaration
template <class t_index>
void construct(t_index &idx, const std::string &file, sdsl::cache_config &config,
               uint8_t num_bytes, trie_tag);

template<trie_construct_algo_type>
class trie_construct {
public:
	typedef sdsl::int_vector<>::size_type           size_type;
	typedef sdsl::int_vector<8>                     text_type;
	typedef sdsl::int_vector_buffer<8>              text_buf_type;
	typedef sdsl::csa_wt<sdsl::wt_blcd<>,0xFFFFFFFF,0xFFFFFFFF> fm_index_type;

public:
	//! prepares input and constructs an fm index. returns length of longest line (lmax)
	static size_type prepare_input( const std::string &file, sdsl::cache_config &config, uint8_t num_bytes, fm_index_type &fm_index );

	//! main construction method, requires text, bwt and suffix array to reside in cache
	//! lmax: length of longest string
	template<class trie_t>
	static void construct( sdsl::cache_config &config, trie_t &trie, fm_index_type &&fm_index, size_type lmax );

private:
	//OHLEBUSCH algorithms
	
	//! constructs mr and cntc
	template<class cntc_t>
	static void construct_mr_cntc( const fm_index_type &fm_index, size_type m, sdsl::bit_vector &mr, cntc_t &cntc );

	//! constructs failure links in form of the balanced parentheses sequence P
	template<class cntc_t>
	static void construct_P( sdsl::bit_vector &mr, cntc_t &cntc, sdsl::bit_vector &P );

	//! constructs the shape (LX and dout): replaces BWT in cache and replaces mr by dout
	static void construct_shape( const fm_index_type &fm_index, sdsl::cache_config &config, sdsl::bit_vector &mr );

	//OHLEBUSCH lightweight algorithms

	//! constructs mr
	static void construct_mr( const fm_index_type &fm_index, size_type m, sdsl::bit_vector &mr );

	//! constructs P in a lightweight manner
	template<class mr_rank_type, class cntc_t>
	static void construct_P_lw( const fm_index_type &fm_index, size_type m, mr_rank_type &mr_rank, cntc_t &cntc, sdsl::bit_vector &P );

	//TXBWT ALGORITHMS
	template<class cntc_t>
	static void construct_mr_cntc_dout( const fm_index_type &fm_index, size_type m, sdsl::bit_vector &mr, cntc_t &cntc, sdsl::bit_vector &dout );

	//! constructs failure links in form of the balanced parentheses sequence P
	//! regarding the prefix interval markings
	template<class cntc_t>
	static void construct_P_txbwt( sdsl::bit_vector &mr, cntc_t &cntc, sdsl::bit_vector &dout, sdsl::bit_vector &P );

	//! constructs the shape (LX, dout and din): replaces BWT in cache and replaces mr by dout
	static void construct_shape_txbwt( const fm_index_type &fm_index, sdsl::cache_config &config, sdsl::bit_vector &mr, sdsl::bit_vector &dout, sdsl::bit_vector &din );

	//MISCELLANOUS ALGORITHMS

	//! constructs the record array R for m leaves of the trie
	static void construct_R( fm_index_type &&fm_index, size_type m, sdsl::cache_config &config, sdsl::int_vector<> &R );

	//! initializes the final extended BWT
	template<class xbwt_t>
	static void init_xbwt( xbwt_t &trie, size_type m, size_type n, sdsl::cache_config &config, sdsl::bit_vector &&dout, sdsl::bit_vector &&P, sdsl::int_vector<> &&R );

	//! initializes the final tunneled extended BWT
	template<class txbwt_t>
	static void init_txbwt( txbwt_t &trie, size_type m, size_type n, sdsl::cache_config &config, sdsl::bit_vector &&dout, sdsl::bit_vector &&din, sdsl::bit_vector &&P, sdsl::int_vector<> &&R );
};

template <typename... Args>
void construct(trie_xbwt<Args...> &idx, const std::string &file, sdsl::cache_config &config,
               uint8_t num_bytes, trie_tag) {
	//prepare input
	typename trie_construct<XBWT>::fm_index_type fm_index;
	auto lmax = trie_construct<XBWT>::prepare_input( file, config, num_bytes, fm_index );

	//construct tries
	{
		auto event = sdsl::memory_monitor::event("TRIECONSTRUCT");
		switch (trie_construct_config::algo) {
		case XBWT:	trie_construct<XBWT>::construct( config, idx, std::move(fm_index), lmax );
				break;
		case XBWT_SC:	trie_construct<XBWT_SC>::construct( config, idx, std::move(fm_index), lmax );
				break;
		case XBWT_LW:	trie_construct<XBWT_LW>::construct( config, idx, std::move(fm_index), lmax );
				break;
		case XBWT_LW_SC:trie_construct<XBWT_LW_SC>::construct( config, idx, std::move(fm_index), lmax );
				break;
		default:	throw std::domain_error("unknown xbwt construction algorithm");
				break;
		}
	}
	//clean up
	if (config.delete_files) {
		sdsl::util::delete_all_files(config.file_map);
	}
};

template <typename... Args>
void construct(trie_txbwt<Args...> &idx, const std::string &file, sdsl::cache_config &config,
               uint8_t num_bytes, trie_tag) {
	//prepare input
	typename trie_construct<TXBWT>::fm_index_type fm_index;
	auto lmax = trie_construct<TXBWT>::prepare_input( file, config, num_bytes, fm_index );

	//construct tries
	{
		auto event = sdsl::memory_monitor::event("TRIECONSTRUCT");
		switch (trie_construct_config::algo) {
		case TXBWT:	trie_construct<TXBWT>::construct( config, idx, std::move(fm_index), lmax );
				break;
		case TXBWT_SC:	trie_construct<TXBWT_SC>::construct( config, idx, std::move(fm_index), lmax );
				break;
		default:	throw std::domain_error("unknown txbwt construction algorithm");
				break;
		}
	}
	//clean up
	if (config.delete_files) {
		sdsl::util::delete_all_files(config.file_map);
	}
};

//// PREPARATION //////////////////////////////////////////////////////////////

template<trie_construct_algo_type a>
typename trie_construct<a>::size_type trie_construct<a>::prepare_input( const std::string &file,
		sdsl::cache_config &config, uint8_t num_bytes, fm_index_type &fm_index ) {
	if (num_bytes != 1) {
		throw std::invalid_argument("tries can be constructed only from plain byte files");
	}
	const char *KEY_TEXT = sdsl::key_text_trait<8>::KEY_TEXT;

	//// prepare input text ///////////////////////////////////////////////
	size_type lmax = 0;
	{
		auto event = sdsl::memory_monitor::event("parse input text");
		text_type text;
		load_vector_from_file(text, file, 1);
		if (text[text.size() - 1] != '\n') { //append an empty line at end
			text.resize(text.size()+1);
			text[text.size()-1] = '\n';
		}
		//reverse lines and replace newlines by nullbytes
		size_type linestart = 0;
		for (size_type i = 0; i < text.size(); i++) {
			if (text[i] == '\n') {
				if (i - linestart > lmax) {
					lmax = i - linestart;
				}

				text[i] = '\0';
				size_type lineend = i;
				while (linestart < lineend--) {
					auto tmp = text[lineend];
					text[lineend] = text[linestart];
					text[linestart] = tmp;
					++linestart;
				}
				linestart = i+1;

			} else if (text[i] == '\0') {
				throw std::invalid_argument("Input contains nullbytes");
			}
		}

		sdsl::store_to_cache(text, KEY_TEXT, config);
		sdsl::register_cache_file(KEY_TEXT, config);
	}
	//// construct FM index ///////////////////////////////////////////////
	{
		auto event = sdsl::memory_monitor::event("FMINDEX");
		bool delete_files = config.delete_files;
		config.delete_files = false; //need suffix array later
		sdsl::construct( fm_index, file, config, 1 );
		config.delete_files = delete_files;
	}
	return lmax;
};

//// XBWT INITIALIZATION //////////////////////////////////////////////////////
template<trie_construct_algo_type a>
template<class xbwt_t>
void trie_construct<a>::init_xbwt( xbwt_t &trie, size_type m, size_type n, sdsl::cache_config &config, sdsl::bit_vector &&dout, sdsl::bit_vector &&P, sdsl::int_vector<> &&R ) {
	trie.n_leaves = m;
	trie.n_strings = n;

	//set up wavelet tree
	text_buf_type LX( cache_file_name(sdsl::key_bwt_trait<8>::KEY_BWT, config) );
	typename xbwt_t::wt_type LX_wt( LX, dout.size() - 1 );
	std::swap( LX_wt, trie.L );
	
	//set up C array
	trie.C.resize( 256 );
	trie.C[0] = 0;
	size_type sum = 1; //count only 1 nullbyte
	for (size_type c = 1; c < 256u; c++) {
		trie.C[c] = sum;
		sum += trie.L.rank( trie.L.size(), c );
	}

	// construct rank and select support for Dout
	typename xbwt_t::bit_vector_type dout_bv( std::move(dout) );
	std::swap( dout_bv, trie.Dout );
	sdsl::util::init_support(trie.Dout_rank, &trie.Dout);
	sdsl::util::init_support(trie.Dout_select, &trie.Dout);

	//construct support for P
	std::swap( P, trie.P );
	sdsl::util::init_support(trie.P_bps_support, &trie.P);
	
	//copy R array to trie
	std::swap( R, trie.R );
};

//// TXBWT INITIALIZATION /////////////////////////////////////////////////////
template<trie_construct_algo_type a>
template<class txbwt_t>
void trie_construct<a>::init_txbwt( txbwt_t &trie, size_type m, size_type n, sdsl::cache_config &config, sdsl::bit_vector &&dout, sdsl::bit_vector &&din, sdsl::bit_vector &&P, sdsl::int_vector<> &&R ) {
	//do normal xbwt initialization
	init_xbwt( trie, m, n, config, std::move(dout), std::move(P), std::move(R) );

	// construct rank and select support for Din
	typename txbwt_t::bit_vector_type din_bv( std::move(din) );
	std::swap( din_bv, trie.Din );
	sdsl::util::init_support(trie.Din_rank, &trie.Din);
	sdsl::util::init_support(trie.Din_select, &trie.Din);
};
	

//// OHLEBUSCH ////////////////////////////////////////////////////////////////

template<>
template<class trie_t>
void trie_construct<XBWT>::construct( sdsl::cache_config &config, trie_t &trie, fm_index_type &&fm_index,
                                     size_type lmax ) {
	size_type m = fm_index.wavelet_tree.rank( fm_index.size(), '\0');
	size_type n = fm_index.size();
	sdsl::bit_vector mr( fm_index.size() + 1, 0 );
	sdsl::bit_vector P;
	sdsl::int_vector<> R;

	{
		sdsl::int_vector<> cntc(fm_index.size(), 0, sdsl::bits::hi(lmax+1) + 1);
		construct_mr_cntc( fm_index, m, mr, cntc );
		construct_P( mr, cntc, P );
	}
	construct_shape( fm_index, config,  mr );
	construct_R( std::move(fm_index), m, config, R );
	init_xbwt( trie, m, n, config, std::move(mr), std::move(P), std::move(R) );
};

template<trie_construct_algo_type a>
template<class cntc_t>
void trie_construct<a>::construct_mr_cntc( const fm_index_type &fm_index, size_type m, sdsl::bit_vector &mr, cntc_t &cntc ) {
	auto &wt = fm_index.wavelet_tree;

	//variables for interval symbols
	typename fm_index_type::wavelet_tree_type::size_type k;
	std::vector<typename fm_index_type::wavelet_tree_type::value_type> cs(fm_index.sigma);
	std::vector<typename fm_index_type::wavelet_tree_type::size_type> rank_c_lb(fm_index.sigma);
	std::vector<typename fm_index_type::wavelet_tree_type::size_type> rank_c_rb(fm_index.sigma);

	//initialization
	mr[wt.size()] = 1;
	std::queue<std::tuple<size_type, size_type, size_type>> Q;
	Q.emplace( 0u, m, wt.size() );

	while (!Q.empty()) {
		size_type i_d = std::get<0>(Q.front());
		size_type j_d = std::get<1>(Q.front());
		size_type j   = std::get<2>(Q.front());
		Q.pop();

		mr[i_d] = 1;
		cntc[j-1]++;

		sdsl::interval_symbols(wt, i_d, j_d, k, cs, rank_c_lb, rank_c_rb);
		for (size_type l = 0; l < k; l++) {
			auto c = cs[l];
			if (c != '\0') {
				size_type lb_d = fm_index.C[fm_index.char2comp[c]] + rank_c_lb[l];
				size_type rb_d = fm_index.C[fm_index.char2comp[c]] + rank_c_rb[l];
				size_type rb = (j_d == j) 
				             ? rb_d
				             : fm_index.C[fm_index.char2comp[c]] + wt.rank( j, c );
				
				Q.emplace( lb_d, rb_d, rb );
			}
		}
	}
};

template<trie_construct_algo_type a>
template<class cntc_t>
void trie_construct<a>::construct_P( sdsl::bit_vector &mr, cntc_t &cntc, sdsl::bit_vector &P ) {
	//get number of nodes
	size_type n_t = 0;
	for (size_type i = 0; i < cntc.size(); i++) {
		n_t+=mr[i];
	}
	P.resize( 2 * n_t );

	//construct sequence
	size_type k = 0;
	for (size_type i = 0; i < cntc.size(); i++) {
		if (mr[i]) {
			P[k++] = 1;
		}
		size_type cc = cntc[i];
		for (size_type j = 0; j < cc; j++) {
			P[k++] = 0;
		}
	}
};

template<trie_construct_algo_type a>
void trie_construct<a>::construct_shape( const fm_index_type &fm_index, sdsl::cache_config &config, sdsl::bit_vector &mr ) {
	auto &wt = fm_index.wavelet_tree;

	//variables for interval symbols
	typename fm_index_type::wavelet_tree_type::size_type m;
	std::vector<typename fm_index_type::wavelet_tree_type::value_type> cs(fm_index.sigma);
	std::vector<typename fm_index_type::wavelet_tree_type::size_type> rank_c_lb(fm_index.sigma);
	std::vector<typename fm_index_type::wavelet_tree_type::size_type> rank_c_rb(fm_index.sigma);

	//overwrite BWT in cache
	text_buf_type LX( cache_file_name(sdsl::key_bwt_trait<8>::KEY_BWT, config) );
	LX.reset();

	size_type i = 0;
	size_type k = 0;
	for (size_type j = 1; j < mr.size(); j++) {
		if (mr[j]) {
			sdsl::interval_symbols(wt, i, j, m, cs, rank_c_lb, rank_c_rb);
			for (size_type l = 0; l < m; l++) {
				LX.push_back( cs[l] );
				mr[k++] = 0;
			}
			mr[k-m] = 1;
			i = j;
		}
	}
	mr[k++] = 1;
	mr.resize( k );
};

template<trie_construct_algo_type a>
void trie_construct<a>::construct_R( fm_index_type &&fm_index, size_type m, sdsl::cache_config &config, sdsl::int_vector<> &R ) {
	//clear R at first moment
	R.width( sdsl::bits::hi( m+1 ) + 1 );
	R.resize( 0 );

	size_type bwt_idx = fm_index.isa[0];
	size_type cyclic_d = fm_index.wavelet_tree.rank( bwt_idx, '\0' ); //compute which $ in L corresponds to first $ in F
	fm_index = fm_index_type(); //get rid of fm index

	//set up an array containing all information
	std::vector<std::pair<size_type,size_type>> A( m );
	{
		sdsl::int_vector_buffer<> SA( cache_file_name(sdsl::conf::KEY_SA, config) );
		size_type j = 1; //sa access position for all except of the cyclic dollar
		for (size_type i = 0; i < m; i++) {
			if (i != cyclic_d) { //assert suffix array position to each $ in L
				A[i] = std::make_pair( SA[j], i );
				++j;
			} else { //assert cyclic position in the special case
				A[i] = std::make_pair( 0, i );
			}
		}
	}

	//sort array and keep only the second key
	std::sort( A.begin(), A.end() );
	for (size_type i = 0; i < m; i++) {
		if (i % 2 == 0) {
			A[i/2].first = A[i].second;
		} else {
			A[i/2].second = A[i].second;
		}
	}
	A.resize( m / 2 + 1 ); //make place for R array

	//copy values into R
	R.resize( m+1 );
	for (size_type i = 0; i < m; i++) {
		size_type j = (i % 2 == 0) ? A[i/2].first : A[i/2].second;
		R[j] = i;
	}
	R[m] = m; //dummy entry	
};

//// OHLEBUSCH SUCCINCT COUNTER ///////////////////////////////////////////////

template<>
template<class trie_t>
void trie_construct<XBWT_SC>::construct( sdsl::cache_config &config, trie_t &trie, fm_index_type &&fm_index,
                                     SDSL_UNUSED size_type lmax ) {
	size_type m = fm_index.wavelet_tree.rank( fm_index.size(), '\0');
	size_type n = fm_index.size();
	sdsl::bit_vector mr( fm_index.size() + 1, 0 );
	sdsl::bit_vector P;
	sdsl::int_vector<> R;

	{
		succinct_counter<> cntc(fm_index.size());
		construct_mr_cntc( fm_index, m, mr, cntc );
		construct_P( mr, cntc, P );
	}
	construct_shape( fm_index, config,  mr );
	construct_R( std::move(fm_index), m, config, R );
	init_xbwt( trie, m, n, config, std::move(mr), std::move(P), std::move(R) );
};

//// OHLEBUSCH LIGHTWEIGHT ////////////////////////////////////////////////////

template<>
template<class trie_t>
void trie_construct<XBWT_LW>::construct( sdsl::cache_config &config, trie_t &trie, fm_index_type &&fm_index,
                                     size_type lmax ) {
	size_type m = fm_index.wavelet_tree.rank( fm_index.size(), '\0');
	size_type n = fm_index.size();
	sdsl::bit_vector mr( fm_index.size() + 1, 0 );
	sdsl::bit_vector P;
	sdsl::int_vector<> R;

	{
		construct_mr( fm_index, m, mr );
		sdsl::bit_vector::rank_1_type mr_rank(&mr);
		size_type num_nodes = mr_rank(mr.size() - 1 );
		sdsl::int_vector<> cntc( num_nodes, 0, sdsl::bits::hi(lmax+1) + 1 );
		construct_P_lw( fm_index, m, mr_rank, cntc, P );
	}
	construct_shape( fm_index, config,  mr );
	construct_R( std::move(fm_index), m, config, R );
	init_xbwt( trie, m, n, config, std::move(mr), std::move(P), std::move(R) );
};

template<trie_construct_algo_type a>
void trie_construct<a>::construct_mr( const fm_index_type &fm_index, size_type m, sdsl::bit_vector &mr ) {
	auto &wt = fm_index.wavelet_tree;

	//variables for interval symbols
	typename fm_index_type::wavelet_tree_type::size_type k;
	std::vector<typename fm_index_type::wavelet_tree_type::value_type> cs(fm_index.sigma);
	std::vector<typename fm_index_type::wavelet_tree_type::size_type> rank_c_lb(fm_index.sigma);
	std::vector<typename fm_index_type::wavelet_tree_type::size_type> rank_c_rb(fm_index.sigma);

	//initialization
	mr[wt.size()] = 1;
	std::queue<std::tuple<size_type, size_type, size_type>> Q;
	Q.emplace( 0u, m, wt.size() );

	while (!Q.empty()) {
		size_type i_d = std::get<0>(Q.front());
		size_type j_d = std::get<1>(Q.front());
		size_type j   = std::get<2>(Q.front());
		Q.pop();

		mr[i_d] = 1;

		sdsl::interval_symbols(wt, i_d, j_d, k, cs, rank_c_lb, rank_c_rb);
		for (size_type l = 0; l < k; l++) {
			auto c = cs[l];
			if (c != '\0') {
				size_type lb_d = fm_index.C[fm_index.char2comp[c]] + rank_c_lb[l];
				size_type rb_d = fm_index.C[fm_index.char2comp[c]] + rank_c_rb[l];
				size_type rb = (j_d == j) 
				             ? rb_d
				             : fm_index.C[fm_index.char2comp[c]] + wt.rank( j, c );
				
				Q.emplace( lb_d, rb_d, rb );
			}
		}
	}
};

template<trie_construct_algo_type a>
template<class mr_rank_type, class cntc_t>
void trie_construct<a>::construct_P_lw( const fm_index_type &fm_index, size_type m, mr_rank_type &mr_rank, cntc_t &cntc, sdsl::bit_vector &P ) {
	auto &wt = fm_index.wavelet_tree;

	{
		//variables for interval symbols
		typename fm_index_type::wavelet_tree_type::size_type k;
		std::vector<typename fm_index_type::wavelet_tree_type::value_type> cs(fm_index.sigma);
		std::vector<typename fm_index_type::wavelet_tree_type::size_type> rank_c_lb(fm_index.sigma);
		std::vector<typename fm_index_type::wavelet_tree_type::size_type> rank_c_rb(fm_index.sigma);

		//initialization
		std::queue<std::tuple<size_type, size_type, size_type>> Q;
		Q.emplace( 0u, m, wt.size() );

		while (!Q.empty()) {
			size_type i_d = std::get<0>(Q.front());
			size_type j_d = std::get<1>(Q.front());
			size_type j   = std::get<2>(Q.front());
			Q.pop();

			cntc[mr_rank(j) - 1]++; // count endpositions (not inclusive)

			sdsl::interval_symbols(wt, i_d, j_d, k, cs, rank_c_lb, rank_c_rb);
			for (size_type l = 0; l < k; l++) {
				auto c = cs[l];
				if (c != '\0') {
					size_type lb_d = fm_index.C[fm_index.char2comp[c]] + rank_c_lb[l];
					size_type rb_d = fm_index.C[fm_index.char2comp[c]] + rank_c_rb[l];
					size_type rb = (j_d == j) 
						     ? rb_d
						     : fm_index.C[fm_index.char2comp[c]] + wt.rank( j, c );
				
					Q.emplace( lb_d, rb_d, rb );
				}
			}
		}
	}
	{
		//compute P
		P.resize(2 * cntc.size());
		size_type k = 0;
		for (size_type i = 0; i < cntc.size(); i++) {
			P[k++] = 1;
			size_type cc = cntc[i];
			for (size_type j = 0; j < cc; j++) {
				P[k++] = 0;
			}
		}
	}
};

//// OHLEBUSCH LIGHTWEIGHT SUCCINCT COUNTER ///////////////////////////////////

template<>
template<class trie_t>
void trie_construct<XBWT_LW_SC>::construct( sdsl::cache_config &config, trie_t &trie, fm_index_type &&fm_index,
                                     SDSL_UNUSED size_type lmax ) {
	size_type m = fm_index.wavelet_tree.rank( fm_index.size(), '\0');
	size_type n = fm_index.size();
	sdsl::bit_vector mr( fm_index.size() + 1, 0 );
	sdsl::bit_vector P;
	sdsl::int_vector<> R;

	{
		construct_mr( fm_index, m, mr );
		sdsl::bit_vector::rank_1_type mr_rank(&mr);
		size_type num_nodes = mr_rank(mr.size() - 1 );
		succinct_counter<> cntc( num_nodes );
		construct_P_lw( fm_index, m, mr_rank, cntc, P );
	}
	construct_shape( fm_index, config,  mr );
	construct_R( std::move(fm_index), m, config, R );
	init_xbwt( trie, m, n, config, std::move(mr), std::move(P), std::move(R) );
};

//// TXBWT ////////////////////////////////////////////////////////////////////

template<>
template<class trie_t>
void trie_construct<TXBWT>::construct( sdsl::cache_config &config, trie_t &trie, fm_index_type &&fm_index,
                                     size_type lmax ) {
	size_type m = fm_index.wavelet_tree.rank( fm_index.size(), '\0');
	size_type n = fm_index.size();
	sdsl::bit_vector mr( fm_index.size() + 1, 0 );
	sdsl::bit_vector dout( fm_index.size() + 1, 0 );
	sdsl::bit_vector din( fm_index.size() + 1, 0 );
	sdsl::bit_vector P;
	sdsl::int_vector<> R;

	{
		//find good order for kmer prefix intervals to be tunneled
		dbg_algorithms::find_min_dbg( fm_index, dout, config );
		din = dout;

		//construct XBWT construction components and set additional 1-bits
		//in dout where a node has outdegree greater one (omits conflicts between tunnels
		//and branching nodes
		sdsl::int_vector<> cntc(fm_index.size(), 0, sdsl::bits::hi(lmax+1) + 1);
		construct_mr_cntc_dout( fm_index, m, mr, cntc, dout );

		//set up a marking of tunnelable non-conflicting k-mer prefix intervals
		for (size_type i = 0; i < m; i++) {
			dout[i] = 1; //avoid removal of $'s
		}
		dbg_algorithms::mark_prefix_intervals( fm_index, dout, din );
		
		//construct failure links
		construct_P_txbwt( mr, cntc, dout, P );
	}
	//and shape
	construct_shape_txbwt( fm_index, config, mr, dout, din );
	construct_R( std::move(fm_index), m, config, R );	
	init_txbwt( trie, m, n, config, std::move(dout), std::move(din), std::move(P), std::move(R) );
};

template<trie_construct_algo_type a>
template<class cntc_t>
void trie_construct<a>::construct_mr_cntc_dout( const fm_index_type &fm_index, size_type m, sdsl::bit_vector &mr, cntc_t &cntc, sdsl::bit_vector &dout ) {
	auto &wt = fm_index.wavelet_tree;

	//variables for interval symbols
	typename fm_index_type::wavelet_tree_type::size_type k;
	std::vector<typename fm_index_type::wavelet_tree_type::value_type> cs(fm_index.sigma);
	std::vector<typename fm_index_type::wavelet_tree_type::size_type> rank_c_lb(fm_index.sigma);
	std::vector<typename fm_index_type::wavelet_tree_type::size_type> rank_c_rb(fm_index.sigma);

	//initialization
	mr[wt.size()] = 1;
	std::queue<std::tuple<size_type, size_type, size_type>> Q;
	Q.emplace( 0u, m, wt.size() );

	while (!Q.empty()) {
		size_type i_d = std::get<0>(Q.front());
		size_type j_d = std::get<1>(Q.front());
		size_type j   = std::get<2>(Q.front());
		Q.pop();

		mr[i_d] = 1;
		cntc[j-1]++;

		sdsl::interval_symbols(wt, i_d, j_d, k, cs, rank_c_lb, rank_c_rb);
		if (k > 1u)	dout[i_d+1] = 1; //mark nodes with more than one outgoing edge
		for (size_type l = 0; l < k; l++) {
			auto c = cs[l];
			if (c != '\0') {
				size_type lb_d = fm_index.C[fm_index.char2comp[c]] + rank_c_lb[l];
				size_type rb_d = fm_index.C[fm_index.char2comp[c]] + rank_c_rb[l];
				size_type rb = (j_d == j) 
				             ? rb_d
				             : fm_index.C[fm_index.char2comp[c]] + wt.rank( j, c );
				
				Q.emplace( lb_d, rb_d, rb );
			}
		}
	}
};

template<trie_construct_algo_type a>
template<class cntc_t>
void trie_construct<a>::construct_P_txbwt( sdsl::bit_vector &mr, cntc_t &cntc, sdsl::bit_vector &dout, sdsl::bit_vector &P ) {

	//get number of failure links
	size_type n_fl = 0;
	for (size_type i = 0; i < cntc.size(); i++) {
		n_fl += (mr[i] & dout[i]);
	}
	P.resize( 2 * n_fl );

	//construct sequence
	size_type k = 0;
	size_type skip = 0;
	for (size_type i = 0; i < cntc.size(); i++) {
		//write opening parenthesis in case of explicit or implicit failure link
		if (mr[i]) {
			if (dout[i] == 1) {
				P[k++] = 1;
			} else {
				++skip;
			}
		}

		//write closing parentheses
		size_type cc = cntc[i];
		if (cc > skip) {
			cc -= skip;
			for (size_type j = 0; j < cc; j++) {
				P[k++] = 0;
			}
			skip = 0;
		} else {
			skip -= cc;
		}
	}
}

template<trie_construct_algo_type a>
void trie_construct<a>::construct_shape_txbwt( const fm_index_type &fm_index, sdsl::cache_config &config,
			sdsl::bit_vector &mr, sdsl::bit_vector &dout, sdsl::bit_vector &din ) {
	auto &wt = fm_index.wavelet_tree;

	//variables for interval symbols
	typename fm_index_type::wavelet_tree_type::size_type m;
	std::vector<typename fm_index_type::wavelet_tree_type::value_type> cs(fm_index.sigma);
	std::vector<typename fm_index_type::wavelet_tree_type::size_type> rank_c_lb(fm_index.sigma);
	std::vector<typename fm_index_type::wavelet_tree_type::size_type> rank_c_rb(fm_index.sigma);

	//overwrite BWT in cache
	text_buf_type LX( cache_file_name(sdsl::key_bwt_trait<8>::KEY_BWT, config) );
	LX.reset();

	size_type i = 0;
	size_type i_out = 0;
	size_type i_in = 0;
	for (size_type j = 1; j < mr.size(); j++) {
		if (mr[j]) {
			bool b = dout[i];
			if (din[i] == 1) {
				sdsl::interval_symbols(wt, i, j, m, cs, rank_c_lb, rank_c_rb);
				for (size_type l = 0; l < m; l++) {
					LX.push_back( cs[l] );
					dout[i_out++] = 0;
				}
				dout[i_out-m] = b;
			}
			if (b) { //dout[i] = 1
				din[i_in++] = din[i];
			}
			i = j;
		}
	}
	dout[i_out++] = 1;	dout.resize( i_out );
	din[i_in++] = 1;	din.resize( i_in );
}

//// TXBWT SUCCINCT COUNTER ///////////////////////////////////////////////////

template<>
template<class trie_t>
void trie_construct<TXBWT_SC>::construct( sdsl::cache_config &config, trie_t &trie, fm_index_type &&fm_index,
                                     SDSL_UNUSED size_type lmax ) {
	size_type m = fm_index.wavelet_tree.rank( fm_index.size(), '\0');
	size_type n = fm_index.size();
	sdsl::bit_vector mr( fm_index.size() + 1, 0 );
	sdsl::bit_vector dout( fm_index.size() + 1, 0 );
	sdsl::bit_vector din( fm_index.size() + 1, 0 );
	sdsl::bit_vector P;
	sdsl::int_vector<> R;

	{
		//find good order for kmer prefix intervals to be tunneled
		dbg_algorithms::find_min_dbg( fm_index, dout, config );
		din = dout;

		//construct XBWT construction components and set additional 1-bits
		//in dout where a node has outdegree greater one (omits conflicts between tunnels
		//and branching nodes
		succinct_counter<> cntc(fm_index.size());
		construct_mr_cntc_dout( fm_index, m, mr, cntc, dout );

		//set up a marking of tunnelable non-conflicting k-mer prefix intervals
		for (size_type i = 0; i < m; i++) {
			dout[i] = 1; //avoid removal of $'s
		}
		dbg_algorithms::mark_prefix_intervals( fm_index, dout, din );
		
		//construct failure links
		construct_P_txbwt( mr, cntc, dout, P );
	}
	//and shape
	construct_shape_txbwt( fm_index, config, mr, dout, din );
	construct_R( std::move(fm_index), m, config, R );	
	init_txbwt( trie, m, n, config, std::move(dout), std::move(din), std::move(P), std::move(R) );
};

#endif

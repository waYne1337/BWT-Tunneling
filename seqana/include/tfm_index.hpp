/*
 * tfm_index.hpp for BWT Tunneling
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
 * tfm_index.hpp for Edge minimization in de Bruijn graphs
 * Copyright (c) 2019 Uwe Baier, Pascal Weber All Rights Reserved.
 */

#ifndef TFM_INDEX_HPP
#define TFM_INDEX_HPP

#include <assert.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>
#include <sdsl/construct_lcp_helper.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/select_support.hpp>
#include <sdsl/util.hpp>
#include <sdsl/wavelet_trees.hpp>

#include <algorithm>
#include <limits>
#include <utility>

#include "dbg_algorithms.hpp"

struct tfm_index_tag {};

//! a class representing a tunneled fm-index
template<class t_wt_type =       sdsl::wt_blcd<>,
         class t_bv_type =       typename t_wt_type::bit_vector_type,
         class t_rank_type =     typename t_bv_type::rank_1_type,
         class t_select_type =   typename t_bv_type::select_1_type>
class tfm_index {
public:
	typedef tfm_index_tag                           index_category;
	typedef sdsl::byte_alphabet_tag                 alphabet_category;
	typedef sdsl::int_vector<>::size_type           size_type;

	typedef sdsl::int_vector<8>                     text_type;
	typedef typename t_wt_type::value_type          value_type;

	typedef t_wt_type                               wt_type;
        typedef t_bv_type                               bit_vector_type;
        typedef t_rank_type                             rank_type;
	typedef t_select_type                           select_type;

	//first index is next outgoing edge, second index is tunnel entry offset
	typedef std::pair<size_type,size_type>          nav_type;

private:
	template <typename t_tfm_index_type>
	friend void construct_tfm_index( t_tfm_index_type &tfm_index, uint64_t text_len, 
		sdsl::int_vector_buffer<8> &&L_buf, sdsl::bit_vector &&dout, sdsl::bit_vector &&din );

	size_type                                       text_len; //original textlen
	wt_type                                         m_L;
	std::vector<size_type>                          m_C;
	bit_vector_type                                 m_dout;
	rank_type                                       m_dout_rank;
	select_type                                     m_dout_select;
        bit_vector_type                                 m_din;
	rank_type                                       m_din_rank;
	select_type                                     m_din_select;

public:
	const wt_type &                                        L = m_L;
	const std::vector<size_type> &                         C = m_C;
        const bit_vector_type &                                dout = m_dout;
	const rank_type &                                      dout_rank = m_dout_rank;
	const select_type &                                    dout_select = m_dout_select;
	const bit_vector_type &                                din = m_din;
	const rank_type &                                      din_rank = m_din_rank;
	const select_type &                                    din_select = m_din_select;

	//! returns the size of the original string
	size_type size() const {
		return text_len;
	};

	//! returns the end, i.e. the position in L where the string ends
	nav_type end() const {
		return std::make_pair( (size_type)0, (size_type)0 );
	}

	//! returns the character preceding the current position
	value_type preceding_char( const nav_type &pos ) const {
		return L[pos.first];
	}

	//! Operation performs an backward step from current position.
	//! function sets posm to the new value and returns the result
	//! of preceding_char( pos ) before the backward step was performed
	value_type backwardstep( nav_type &pos ) const {
		size_type &i = pos.first; //create references into position pair
		size_type &o = pos.second;

		//navigate to next entry
		auto is = L.inverse_select( i );
		auto c = is.second;
		i = C[c] + is.first;

		//check for the start of a tunnel
		auto din_rank_ip1 = din_rank( i + 1 );
		if (din[i] == 0) {
			o = i - din_select( din_rank_ip1 ); //save offset to uppermost entry edge
		}
		//navigate to outedges of current node
		i = dout_select( din_rank_ip1 );

		//check for end of a tunnel
		if (dout[i+1] == 0) {
			i += o; //jump back offset
			o = 0;
		}
		return c;
	};

	//! serializes opbject
	size_type serialize(std::ostream &out, sdsl::structure_tree_node *v,
                          std::string name) const {

		sdsl::structure_tree_node *child =
			sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
		size_type written_bytes = 0;
		written_bytes += sdsl::write_member(text_len, out, child, "text_len");

		written_bytes += m_L.serialize(out, child, "L");
		written_bytes += sdsl::serialize(m_C, out, child, "C");

		written_bytes += m_dout.serialize(out, child, "dout");
		written_bytes += m_dout_rank.serialize(out, child, "dout_rank");
		written_bytes += m_dout_select.serialize(out, child, "dout_select");

		written_bytes += m_din.serialize(out, child, "din");
		written_bytes += m_din_rank.serialize(out, child, "din_rank");
		written_bytes += m_din_select.serialize(out, child, "din_select");

		sdsl::structure_tree::add_size(child, written_bytes);
		return written_bytes;
	};

	//! loads a serialized object
	void load(std::istream &in) {

		sdsl::read_member( text_len, in );

		m_L.load(in);
		sdsl::load(m_C, in);

		m_dout.load(in);
		m_dout_rank.load(in, &m_dout);
		m_dout_select.load(in, &m_dout);

		m_din.load(in);
		m_din_rank.load(in, &m_din);
		m_din_select.load(in, &m_din);
	};
};

//// SPECIAL CONSTRUCTION FOR TUNNELED FM INDEX ///////////////////////////////

template <class t_index>
void construct(t_index &idx, const std::string &file, sdsl::cache_config &config,
               uint8_t num_bytes, tfm_index_tag) {
	assert( num_bytes == 1 ); //only byte input is allowed

	//create a normal fm index
	sdsl::csa_wt<sdsl::wt_blcd<>,0xFFFFFFFF,0xFFFFFFFF> csa;
	{

		construct( csa, file, config, 1 );
	}

	//run construction algorithm
	{
		auto event = sdsl::memory_monitor::event("construct tunneled fm index");
		construct_tfm_index( idx, std::move( csa ), config );
	}
};

//! function constructs a tfm index using a compressed suffix array in form of a BWT in a wavelet tree.
//! note that the csa is erased during construction
//! function returns the result of the dbg_algorithms::find_min_dbg - function
template<class t_tfm_index_type,
         class t_csa_wt_type>
std::pair<typename t_tfm_index_type::size_type,typename t_tfm_index_type::size_type>
construct_tfm_index( t_tfm_index_type &tfm_index, t_csa_wt_type &&csa, sdsl::cache_config &config ) {
	typedef typename t_tfm_index_type::size_type size_type;
	std::pair<size_type,size_type> dbg_res;

	//find minimal edge-reduced DBG and store kmer bounds in a bitvector B
        sdsl::bit_vector B;
	{
		auto event = sdsl::memory_monitor::event("FINDMINDBG");
		dbg_res = dbg_algorithms::find_min_dbg( csa, B, config );
	}

	//use bitvector to determine prefix intervals to be tunneled
	auto event = sdsl::memory_monitor::event("TFMINDEXCONSTRUCT");
	sdsl::bit_vector dout = B;
	sdsl::bit_vector din; std::swap( din, B );
	dbg_algorithms::mark_prefix_intervals( csa, dout, din );

	//create a buffer for newly constructed L
	std::string tmp_key = sdsl::util::to_string(sdsl::util::pid())+"_"+sdsl::util::to_string(sdsl::util::id());
	std::string tmp_file_name = sdsl::cache_file_name(tmp_key, config);
	{
		sdsl::int_vector_buffer<8> L_buf(tmp_file_name, std::ios::out);

		//remove redundant entries from L, dout and din
		size_type p = 0;
		size_type q = 0;
		for (size_type i = 0; i < csa.size(); i++) {
			if (din[i] == 1) {
				L_buf.push_back( csa.wavelet_tree[i] );
				dout[p++] = dout[i];
			}
			if (dout[i] == 1) {
				din[q++] = din[i];
			}
		}
		dout[p++] = 1;	din[q++] = 1;
		dout.resize( p );
		din.resize( q );

		uint64_t text_len = csa.size();
		csa = t_csa_wt_type(); //remove csa object as it is no longer required

		construct_tfm_index( tfm_index, text_len, std::move( L_buf ), std::move( dout ), std::move( din ) );
	}
	//remove buffer for L
	sdsl::remove(tmp_file_name);
	return dbg_res;
};

template<class t_tfm_index_type>
void construct_tfm_index( t_tfm_index_type &tfm_index, uint64_t text_len, sdsl::int_vector_buffer<8> &&L_buf, sdsl::bit_vector &&dout, sdsl::bit_vector &&din ) {
	//set original string size
	tfm_index.text_len = text_len;

	//construct tfm index from L, din and dout
	typedef typename t_tfm_index_type::wt_type         wt_type;
	typedef typename t_tfm_index_type::bit_vector_type bv_type;

	//wavelet tree of L
	tfm_index.m_L = wt_type( L_buf, L_buf.size() );
	sdsl::create_C_array( tfm_index.m_C, tfm_index.m_L );

	//dout
	tfm_index.m_dout = bv_type( std::move( dout ) );
	sdsl::util::init_support( tfm_index.m_dout_rank, &tfm_index.m_dout );
	sdsl::util::init_support( tfm_index.m_dout_select, &tfm_index.m_dout );

	//din
	tfm_index.m_din = bv_type( std::move( din ) );
	sdsl::util::init_support( tfm_index.m_din_rank, &tfm_index.m_din );
	sdsl::util::init_support( tfm_index.m_din_select, &tfm_index.m_din );
};

#endif

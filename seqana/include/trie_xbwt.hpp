/*
 * trie_xbwt.hpp for BWT Tunneling
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
 * xbwt.hpp for trickier xbwt tricks
 * Copyright (c) 2018 Stefan Stauß, Uwe Baier All Rights Reserved.
 */

#ifndef TRIE_XBWT_HPP
#define TRIE_XBWT_HPP

#include <algorithm>
#include <utility>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/bits.hpp>
#include <sdsl/bp_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/select_support.hpp>
#include <sdsl/util.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "tries.hpp"

//! xbwt main class
template<class t_wt_type =       sdsl::wt_blcd<>,
         class t_bv_type =       typename t_wt_type::bit_vector_type,
         class t_rank_type =     typename t_bv_type::rank_1_type,
         class t_select_type =   typename t_bv_type::select_1_type,
         class t_bp_support =    sdsl::bp_support_sada<>>
class trie_xbwt {
public:
	typedef trie_tag                                index_category;
	typedef sdsl::byte_alphabet_tag                 alphabet_category;
	typedef sdsl::int_vector<>::size_type           size_type;
	typedef unsigned char                           value_type;
	typedef size_type                               node_type;

	typedef sdsl::int_vector<8>                     text_type;

	typedef t_wt_type                               wt_type;
        typedef t_bv_type                               bit_vector_type;
        typedef t_rank_type                             rank_type;
	typedef t_select_type                           select_type;
	typedef t_bp_support                            bp_support;
private:
	template<trie_construct_algo_type> friend class trie_construct;
	size_type                                       n_leaves;
	size_type                                       n_strings;

	wt_type                                         L;

	bit_vector_type                                 Dout;
	rank_type                                       Dout_rank;
	select_type                                     Dout_select;

	sdsl::int_vector<64>                            C;

	sdsl::bit_vector                                P;
	bp_support                                      P_bps_support;

	sdsl::int_vector<>                              R; //R has size n_leaves+1, last entry in R contains n_leaves
public:
	trie_xbwt() : n_leaves{ 0 } {};

	//! returns number of leaves
	size_type num_leaves() const {
		return n_leaves;
	}

	//! returns the sum of the lengths of all original strings
	size_type strings_length() const {
		return n_strings;
	}

	//! returns id of root node
	node_type root() const {
		return 0;
	};

	//! returns a range of edges for node i (end is exclusive).
	//! number of childs can be computed using second - first
	std::pair<size_type,size_type> get_edge_range(node_type i) const {
		return std::make_pair( Dout_select(i + 1), Dout_select(i + 2) );
	};

	//! returns an empty edge range
	std::pair<size_type,size_type> get_empty_edge_range() const {
		return std::make_pair( (size_type)0, (size_type)0 );
	};

	//! returns a the edge indicated by index e.
	//! function returns a pair of the edge label and an edge identifier.
	//! Combine get_child_range and this function to enumerate all edge labels
	std::pair<value_type,size_type> get_edge(size_type e) const {
		auto edge = L.inverse_select(e);
		return std::make_pair( (value_type)edge.second, (size_type)edge.first );
	};

	//! checks if a edge range contains an outgoing edge labeled c.
	//! function returns a pair of the edge label and an edge identifier in case of success,
	//! a pair containing a nullbyte as character and num_leaves() as identifier otherwise.
	std::pair<value_type,size_type> find_edge( const std::pair<size_type,size_type> &e_range, value_type c ) const {
		size_type r_lb = L.rank(e_range.first, c);
		size_type r_rb = L.rank(e_range.second, c);
		return (r_lb == r_rb)
			? std::make_pair( (value_type)'\0', num_leaves() )
			: std::make_pair( c, r_lb );
	};
	
	//! function returns the record number (line number) of the leaf if edge range
	//! contains an edge to a leaf, or num_leaves() otherwise.
	size_type has_leaf( const std::pair<size_type,size_type> &e_range ) const {
		auto e_id = find_edge( e_range, '\0' );
		return R[e_id.second];
	};
		
	//! computes the node id of the node at which the given edge of outgoing node i points.
	//! Do not use this in case that edge points to a leaf, undefined behaviour!
	node_type follow_edge( const std::pair<value_type,size_type> edge, SDSL_UNUSED node_type i ) const {
		//i is for the interface of a tunneled XBWT here.
		return C[edge.first] + edge.second;
	};

	//! returns the node to which the failure link from node i points to
	//! i mustn't be the root node
	node_type failure_link(node_type i) const {
		size_type k = P_bps_support.select(i + 1);
		size_type j = P_bps_support.enclose(k);
		size_type l = P_bps_support.rank(j);

		return l - 1;
	};

	//TODO: add more functionality if wanted

	//! serializes opbject
	size_type serialize(std::ostream &out, sdsl::structure_tree_node *v,
                          std::string name) const {

		sdsl::structure_tree_node *child =
			sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
		size_type written_bytes = 0;
		written_bytes += sdsl::write_member(n_leaves,out,child, "n_leaves");
		written_bytes += sdsl::write_member(n_strings,out,child, "n_strings");
		written_bytes += L.serialize(out, child, "L");
		written_bytes += Dout.serialize(out, child, "Dout");
		written_bytes += Dout_rank.serialize(out, child, "Dout_rank");
		written_bytes += Dout_select.serialize(out, child, "Dout_select");
		written_bytes += C.serialize(out, child, "C");
		written_bytes += P.serialize(out, child, "P");
		written_bytes += P_bps_support.serialize(out, child, "P_bps_support");
		written_bytes += R.serialize(out, child, "R");
		sdsl::structure_tree::add_size(child, written_bytes);
		return written_bytes;
	};

	//! loads a serialized object
	void load(std::istream &in) {
		sdsl::read_member(n_leaves, in);
		sdsl::read_member(n_strings, in);
		L.load(in);
		Dout.load(in);
		Dout_rank.load(in, &Dout);
		Dout_select.load(in, &Dout);
		C.load(in);
		P.load(in);
		P_bps_support.load(in, &P);
		R.load(in);
	};
};

#endif

/*
 * tries.hpp for BWT Tunneling
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
#ifndef TRIES_HPP
#define TRIES_HPP

#include <stdexcept>

#include <sdsl/construct.hpp>
#include <sdsl/construct_config.hpp>
#include <sdsl/csa_wt.hpp>

//// construction parameters ////
struct trie_tag {};
enum trie_construct_algo_type {XBWT,XBWT_SC,XBWT_LW,XBWT_LW_SC,TXBWT,TXBWT_SC};

class trie_construct_config
{
    public:
        static trie_construct_algo_type algo;

        trie_construct_config() = delete;
};

//set default
trie_construct_algo_type trie_construct_config::algo = XBWT_SC;

#include "trie_xbwt.hpp"
#include "trie_txbwt.hpp"
#include "trie_construct_algorithms.hpp"

//// trie traits ////
template<class trie_t>
struct trie_trait;

template<>
struct trie_trait<trie_xbwt<>> {
	static constexpr const char *name = "TRIE_XBWT<>";
	typedef trie_xbwt<> type;
};

template<>
struct trie_trait<trie_txbwt<>> {
	static constexpr const char *name = "TRIE_TXBWT<>";
	typedef trie_txbwt<> type;
};

#endif

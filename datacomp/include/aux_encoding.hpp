/*
 * aux_encoding.hpp for BWT Tunneling
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
 * aux-encoding.hpp for bwt tunneling
 * Copyright (c) 2017 Uwe Baier All Rights Reserved.
 */

#ifndef AUX_ENCODING_HPP
#define AUX_ENCODING_HPP

#include "bwt_config.hpp"

//! namespace gathering constants for interpretation of the auxiliary data structure
namespace aux_encoding {
	//! regular bwt entry
	const t_size_t REG = 0;
	//! entry indicating the end of a tunnel
	const t_size_t SKP_F = 1;
	//! entry indicating the start of a tunnel
	const t_size_t IGN_L = 2;
	//! entry to be removed
	const t_size_t REM = SKP_F | IGN_L;
	//! alphabet size in auxiliary data structure
	const t_size_t SIGMA = 3;
};

#endif

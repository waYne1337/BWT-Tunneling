/*
 * bcm_poststage.hpp for BWT Tunneling
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
 * bcm-compressor.hpp for bwt tunneling
 * Copyright (c) 2017 Uwe Baier All Rights Reserved.
 */

#ifndef BCM_POSTSTAGE_HPP
#define BCM_POSTSTAGE_HPP

#include "bcm_ss.hpp"

#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>

class bcm_poststage {
	
public:
	//! encodes the transform t
	template<class T>
	static void encode( T &t, std::ostream &out ) {
		bcm::CM cm;
		for (unsigned long i = 0; i < t.size(); i++) {
			cm.Encode( t[i], out );
		}
		cm.Flush(out);
	}

	//! decodes the transform and stores it in t
	template<class T>
	static void decode( std::istream &in, T &t ) {
		bcm::CM cm;
		cm.Init(in);
		for (unsigned long i = 0; i < t.size(); i++) {
			t[i] = cm.Decode(in);
		}
	}
};

#endif

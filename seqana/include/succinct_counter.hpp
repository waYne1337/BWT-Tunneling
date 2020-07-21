/*
 * succinct_counter.hpp for BWT Tunneling
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
//IDEA COMES FROM:
/*
 * Fabio Cunial, Jarno Alanko, and Djamal Belazzougui. “A framework
 * for space-efficient variable-order Markov models.” In: Bioinformatics
 * 35.22 (2019), pp. 4607–4616.
 */

#ifndef SUCCINCT_COUNTER_HPP
#define SUCCINCT_COUNTER_HPP

#include <sdsl/int_vector.hpp>

#include <unordered_map>

//forward declaration
template<uint8_t t_width>
class succinct_counter;

//reference type
template<uint8_t t_width>
class succinct_counter_reference
{
	public:
		typedef typename succinct_counter<t_width>::size_type       size_type;
	private:
		succinct_counter<t_width> &m_counter;
		const size_type m_i;
	public:
		succinct_counter_reference( succinct_counter<t_width> &counter,
		                            const size_type &i ) : m_counter{ counter}, m_i{i} {};

		succinct_counter_reference& operator=(size_type x) {
			m_counter.modify( m_i, x );
			return *this;
		}

		succinct_counter_reference& operator=(const succinct_counter_reference& x) {
			return *this = size_type(x);
		};

		operator size_type() const {
			return m_counter.access(m_i);
		}

	        succinct_counter_reference& operator++() { //prefix
			m_counter.add(m_i, size_type(1));
			return *this;
		}

	        size_type operator++(int) { //postfix
			return m_counter.add(m_i, size_type(1));
		}

	        succinct_counter_reference& operator+=(const size_type x) {
			m_counter.add(m_i, x);
			return *this;
		}

		bool operator==(const succinct_counter_reference& x)const {
			return size_type(*this) == size_type(x);
		}

		bool operator<(const succinct_counter_reference& x)const {
			return size_type(*this) < size_type(x);
		}
};

//class itself
template<uint8_t t_width = 8>
class succinct_counter {
	private:
		static_assert(t_width <= 64 && t_width >= 2, "succinct_counter: width of must be between 2 and 64bits.");
	public:
		typedef sdsl::int_vector<>::size_type       size_type;
        	typedef succinct_counter_reference<t_width> reference;

	private:
		static const size_type max_c_val = (1 << t_width) - 1;

		sdsl::int_vector<t_width> m_c; //counter array for small values
		std::unordered_map<size_type,size_type> m_M; //map for big values

	public:
		succinct_counter(size_type size) : m_c( size ), m_M( size >> t_width ) {}

		bool empty() const	{
			return m_c.empty();
		}

		size_type size() const {
			return m_c.size();
		}

		//low-level functions
		//access entry
		inline size_type access( const size_type &i ) const {
			return (m_c[i] == max_c_val) ? m_M.at(i) : size_type(m_c[i]);
		}

		//set an entry
		inline void modify(const size_type &i, const size_type x ) {
			if (x < max_c_val) {
				if (m_c[i] == max_c_val)	m_M.erase( i );
				m_c[i] = x;
			} else {
				m_c[i] = max_c_val;
				m_M[i] = x;
			}
		}				

		//adds a value and returns value after increment
		inline size_type add( const size_type &i, const size_type x ) {
			auto sum = x + m_c[i];
			if (sum < max_c_val) {
				m_c[i] = sum;
				return sum;
			} else {
				auto &m_val = m_M[i];
				if (m_c[i] < max_c_val) {
					m_c[i] = max_c_val;
					m_val = sum;
				} else {
					m_val += x;
				}
				return m_val;
			}
		}

		//high-level functions with []-operator access
		inline reference operator[](const size_type& i) {
			return reference( *this, i );
		}

		inline size_type operator[](const size_type& i) const {
			return access(i);
		}
};

#endif

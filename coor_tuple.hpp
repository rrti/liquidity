#ifndef LIQUIDITY_COOR_TUPLE_HDR
#define LIQUIDITY_COOR_TUPLE_HDR

#include <algorithm> // std::{min,max}
#include <cassert>
#include <cmath> // std::{sqrt,sin,cos}
#include <cstring> // mem{set,cpy}
#include <initializer_list>

namespace math {
	template<typename type, size_t size>
	struct t_tuple {
	public:
		// note: IL's force {}-syntax in all ctor calls
		t_tuple(const t_tuple& c) { std::memcpy(&a[0], &c.a[0], size * sizeof(type)); }
		t_tuple(const std::initializer_list<type>& args = {}) {
			std::memset(&a[0], 0, sizeof(type) * size);
			std::memcpy(&a[0], args.begin(), std::min(args.size(), size) * sizeof(type));
		}

		t_tuple  operator - () const { t_tuple r; for (size_t n = 0; n < size; n++) r[n] = -(*this)[n]; return r; }

		t_tuple  operator +  (const t_tuple& v) const { t_tuple r; for (size_t n = 0; n < size; n++) r[n] = (*this)[n] + v[n]; return r; }
		t_tuple  operator -  (const t_tuple& v) const { t_tuple r; for (size_t n = 0; n < size; n++) r[n] = (*this)[n] - v[n]; return r; }

		t_tuple& operator += (const t_tuple& v) { for (size_t n = 0; n < size; n++) (*this)[n] += v[n]; return *this; }
		t_tuple& operator -= (const t_tuple& v) { for (size_t n = 0; n < size; n++) (*this)[n] -= v[n]; return *this; }

		t_tuple  operator *  (const type s) const { t_tuple r; for (size_t n = 0; n < size; n++)       r[n]  = (*this)[n] * s; return     r; }
		t_tuple& operator *= (const type s)       {            for (size_t n = 0; n < size; n++) (*this)[n] *=              s; return *this; }

		t_tuple  operator /  (const type s) const { return ((*this) *  (type(1) / s)); }
		t_tuple& operator /= (const type s)       { return ((*this) *= (type(1) / s)); }

		// inner-product
		type operator * (const t_tuple& v) const { type r = type(0); for (size_t n = 0; n < size; n++) r += ((*this)[n] * v[n]); return r; }

		type sq_len() const { return ((*this) * (*this)); }
		type    len() const { return (std::sqrt(sq_len())); }

		type  operator [] (const size_t i) const { assert(i < size); return a[i]; }
		type& operator [] (const size_t i)       { assert(i < size); return a[i]; }

		// convenience accessors
		type  x() const { return (*this)[0]; }
		type  y() const { return (*this)[1]; }
		type  z() const { return (*this)[2]; }
		type  w() const { return (*this)[3]; }
		type& x()       { return (*this)[0]; }
		type& y()       { return (*this)[1]; }
		type& z()       { return (*this)[2]; }
		type& w()       { return (*this)[3]; }

	private:
		type a[size];
	};


	typedef t_tuple<float, 2> t_vec2f;
	typedef t_tuple<float, 3> t_vec3f;
	typedef t_tuple<float, 4> t_vec4f;
};

using math::t_vec3f;
using math::t_vec4f;

#endif


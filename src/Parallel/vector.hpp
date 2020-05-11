//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

#ifndef INCLUDED_SIMPLE_VECTOR_HPP
#define INCLUDED_SIMPLE_VECTOR_HPP

//! \file      vector.hpp
//! \author    B. Dorschner
//! \version   1.1
//! \date      June 5 2019

#include <type_traits>
#include <array>
#include <iostream>

namespace math
{

template<typename T, std::size_t N>
class vector
{
public: // member types

	using data_type = std::array<T,N>;
	using value_type = typename data_type::value_type;
	using size_type = typename data_type::size_type;
	using difference_type = typename data_type::difference_type;
	using reference = typename data_type::reference;
	using const_reference = typename data_type::const_reference;
	using pointer = typename data_type::pointer;
	using const_pointer = typename data_type::const_pointer;
	using iterator = typename data_type::iterator;
	using const_iterator = typename data_type::const_iterator;
	using reverse_iterator = typename data_type::reverse_iterator;
	using const_reverse_iterator = typename data_type::const_reverse_iterator;

public: // ctors

	vector() = default;
	vector(const T& element) { data.fill(element); }
	vector(const vector&) = default;
	vector(vector&&) = default;
	explicit vector(const T* ptr) { for(unsigned int i=0; i<N; ++i) data[i] = ptr[i]; }
	template<typename T2>
	vector(const vector<T2,N>& other) { for(unsigned int i=0; i<N; ++i) data[i] = static_cast<T>(other.data[i]); }
	vector(const std::array<T,N>& _data) : data(_data) {}
	vector(std::array<T,N>&& _data) : data(std::move(_data)) {}
	
	vector& operator=(const vector&) & = default;
	vector& operator=(vector&&) & = default;
	vector& operator=(const T& element) & { data.fill(element); return *this; }

public: // static member functions

	static constexpr std::size_t size() { return N; }
	
public: // member functions

	reference operator[](size_type pos) { return data[pos]; }
	const_reference operator[](size_type pos) const { return data[pos]; }
	
	reference operator()(size_type i) { return data[i]; }
	const_reference operator()(size_type i) const { return data[i]; }
	
	reference front() { return data.front(); }
	const_reference front() const { return data.front(); }
	reference back() { return data.back(); }
	const_reference back() const { return data.back(); }
	reference x() { return data[0]; }
	const_reference x() const { return data[0]; }
	reference y() { return data[1]; }
	const_reference y() const { return data[1]; }
	reference z() { return data[2]; }
	const_reference z() const { return data[2]; }
	reference w() { return data[3]; }
	const_reference w() const { return data[3]; }
	
	
	friend std::ostream& operator<<(std::ostream& os, const vector& v)
	{
		os << "";
		for (unsigned int i=0; i<N; ++i)
		{
			os << v.data[i];
			if (i<N-1) os << " ";
		}
		os << "";
		return os;
	}

public: // iterators

	iterator begin() { return data.begin(); }
	const_iterator begin() const { return data.begin(); }
	const_iterator cbegin() const { return data.cbegin(); }
	iterator end() { return data.end(); }
	const_iterator end() const { return data.end(); }
	const_iterator cend() const { return data.cend(); }
	reverse_iterator rbegin() { return data.rbegin(); }
	const_reverse_iterator rbegin() const { return data.rbegin(); }
	const_reverse_iterator crbegin() const { return data.crbegin(); }
	reverse_iterator rend() { return data.rend(); }
	const_reverse_iterator rend() const { return data.rend(); }
	const_reverse_iterator crend() const { return data.crend(); }
	
public: // arithmetic member functions

	vector& operator+=(const vector& other) { for (unsigned int i=0; i<N; ++i) data[i]+=other[i]; return *this; }
	vector& operator+=(const T& element) { for (unsigned int i=0; i<N; ++i) data[i]+=element; return *this; }
	vector& operator-=(const vector& other) { for (unsigned int i=0; i<N; ++i) data[i]-=other[i]; return *this; }
	vector& operator-=(const T& element) { for (unsigned int i=0; i<N; ++i) data[i]-=element; return *this; }
	vector& operator*=(const vector& other) { for (unsigned int i=0; i<N; ++i) data[i]*=other[i]; return *this; }
	vector& operator*=(const T& element) { for (unsigned int i=0; i<N; ++i) data[i]*=element; return *this; }
	vector& operator/=(const vector& other) { for (unsigned int i=0; i<N; ++i) data[i]/=other[i]; return *this; }
	vector& operator/=(const T& element) { for (unsigned int i=0; i<N; ++i) data[i]/=element; return *this; }

public: // members

	data_type data;
};


template< class T, std::size_t N >
bool operator==( const math::vector<T,N>& lhs, const math::vector<T,N>& rhs ) { return lhs.data == rhs.data; }

template< class T, std::size_t N >
bool operator!=( const math::vector<T,N>& lhs, const math::vector<T,N>& rhs ) { return lhs.data != rhs.data; }

template< class T, std::size_t N >
bool operator<( const math::vector<T,N>& lhs, const math::vector<T,N>& rhs ) { return lhs.data < rhs.data; }

template< class T, std::size_t N >
bool operator<=( const math::vector<T,N>& lhs, const math::vector<T,N>& rhs ) { return lhs.data <= rhs.data; }

template< class T, std::size_t N >
bool operator>( const math::vector<T,N>& lhs, const math::vector<T,N>& rhs ) { return lhs.data > rhs.data; }

template< class T, std::size_t N >
bool operator>=( const math::vector<T,N>& lhs, const math::vector<T,N>& rhs ) { return lhs.data >= rhs.data; }


template<typename T, std::size_t N>
math::vector<T,N> operator-(math::vector<T,N> v) 
{
	for (unsigned int i=0; i<N; ++i) v.data[i] = -v.data[i]; 
	return std::move(v);
}

template<typename U, typename V, std::size_t N>
math::vector<typename std::common_type<U,V>::type,N> operator+(const math::vector<U,N>& lhs, const math::vector<V,N>& rhs) 
{
	math::vector<typename std::common_type<U,V>::type,N> res(lhs);
	for(unsigned int i=0; i<N; ++i) res[i] += rhs[i];
	return std::move(res); 
}

template<typename U, typename V, std::size_t N>
math::vector<typename std::common_type<U,V>::type,N> operator+(const math::vector<U,N>& lhs, const V& rhs) 
{
	math::vector<typename std::common_type<U,V>::type,N> res(lhs);
	for(unsigned int i=0; i<N; ++i) res[i] += rhs;
	return std::move(res); 
}

template<typename U, typename V, std::size_t N>
math::vector<typename std::common_type<U,V>::type,N> operator+(const V& lhs, const math::vector<U,N>& rhs) 
{
	math::vector<typename std::common_type<U,V>::type,N> res(rhs);
	for(unsigned int i=0; i<N; ++i) res[i] += lhs;
	return std::move(res); 
}

template<typename U, typename V, std::size_t N>
math::vector<typename std::common_type<U,V>::type,N> operator-(const math::vector<U,N>& lhs, const math::vector<V,N>& rhs) 
{
	math::vector<typename std::common_type<U,V>::type,N> res(lhs);
	for(unsigned int i=0; i<N; ++i) res[i] -= rhs[i];
	return std::move(res); 
}

template<typename U, typename V, std::size_t N>
math::vector<typename std::common_type<U,V>::type,N> operator-(const math::vector<U,N>& lhs, const V& rhs) 
{
	math::vector<typename std::common_type<U,V>::type,N> res(lhs);
	for(unsigned int i=0; i<N; ++i) res[i] -= rhs;
	return std::move(res); 
}

template<typename U, typename V, std::size_t N>
math::vector<typename std::common_type<U,V>::type,N> operator-(const V& lhs, const math::vector<U,N>& rhs) 
{
	math::vector<typename std::common_type<U,V>::type,N> res(rhs);
	for(unsigned int i=0; i<N; ++i) res[i] -= lhs;
	return std::move(res); 
}

template<typename U, typename V, std::size_t N>
math::vector<typename std::common_type<U,V>::type,N> operator*(const math::vector<U,N>& lhs, const math::vector<V,N>& rhs) 
{
	math::vector<typename std::common_type<U,V>::type,N> res(lhs);
	for(unsigned int i=0; i<N; ++i) res[i] *= rhs[i];
	return std::move(res); 
}

template<typename U, typename V, std::size_t N>
math::vector<typename std::common_type<U,V>::type,N> operator*(const math::vector<U,N>& lhs, const V& rhs) 
{
	math::vector<typename std::common_type<U,V>::type,N> res(lhs);
	for(unsigned int i=0; i<N; ++i) res[i] *= rhs;
	return std::move(res); 
}

template<typename U, typename V, std::size_t N>
math::vector<typename std::common_type<U,V>::type,N> operator*(const V& lhs, const math::vector<U,N>& rhs) 
{
	math::vector<typename std::common_type<U,V>::type,N> res(rhs);
	for(unsigned int i=0; i<N; ++i) res[i] *= lhs;
	return std::move(res); 
}

template<typename U, typename V, std::size_t N>
math::vector<typename std::common_type<U,V>::type,N> operator/(const math::vector<U,N>& lhs, const math::vector<V,N>& rhs) 
{
	math::vector<typename std::common_type<U,V>::type,N> res(lhs);
	for(unsigned int i=0; i<N; ++i) res[i] /= rhs[i];
	return std::move(res); 
}

template<typename U, typename V, std::size_t N>
math::vector<typename std::common_type<U,V>::type,N> operator/(const math::vector<U,N>& lhs, const V& rhs) 
{
	math::vector<typename std::common_type<U,V>::type,N> res(lhs);
	for(unsigned int i=0; i<N; ++i) res[i] /= rhs;
	return std::move(res); 
}

template<typename U, typename V, std::size_t N>
math::vector<typename std::common_type<U,V>::type,N> operator/(const V& lhs, const math::vector<U,N>& rhs) 
{
	math::vector<typename std::common_type<U,V>::type,N> res(rhs);
	for(unsigned int i=0; i<N; ++i) res[i] = lhs/res[i];
	return std::move(res); 
}


template<typename U, typename V, std::size_t N>
typename std::common_type<U,V>::type dot(const math::vector<U,N>& lhs, const math::vector<V,N>& rhs) 
{
	typename std::common_type<U,V>::type res(lhs[0]*rhs[0]); 
	for (unsigned int i=1; i<N; ++i) 
	res+=lhs[i]*rhs[i]; 
	return std::move(res);
}

template<typename U, std::size_t N>
U norm2(const math::vector<U,N>& v) { return dot(v,v); }

template<typename U, std::size_t N>
U norm_inf(const math::vector<U,N>& v) 
{
	using namespace std;
	auto m = std::fabs(v[0]);
	for (unsigned int i=1; i<N; ++i) m = max(m,std::fabs(v[i]));
	return m; 
}


template<typename U, typename V, std::size_t N>
typename std::enable_if<N==3,math::vector<typename std::common_type<U,V>::type,N>>::type cross(const math::vector<U,N>& lhs, const math::vector<V,N>& rhs)
{
	return math::vector<typename std::common_type<U,V>::type,N>({lhs.y()*rhs.z()-lhs.z()*rhs.y(), lhs.z()*rhs.x()-lhs.x()*rhs.z(), lhs.x()*rhs.y()-lhs.y()*rhs.x()});
}

template<typename U, typename V, std::size_t N>
typename std::enable_if<N==2,typename std::common_type<U,V>::type>::type cross(const math::vector<U,N>& lhs, const math::vector<V,N>& rhs)
{
	return lhs.x()*rhs.y()-lhs.y()*rhs.x();
}


} // namespace math


#endif // INCLUDED_SIMPLE_VECTOR_HPP

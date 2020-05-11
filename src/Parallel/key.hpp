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

#ifndef INCLUDED_KEY_HPP
#define INCLUDED_KEY_HPP

//! \file      key.hpp
//! \author    B. Dorschner, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include <array>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "vector.hpp"

namespace decomposition
{


template <int Dim>
struct Key
{
public: // member types

    //using value_type = uint_fast64_t;
    using value_type = unsigned long long int;
    using scalar_coordinate_type=int;

    using float_type = double;
    template<typename U>
    using vector_type = math::vector<U,Dim>;
    using coordinate_type = vector_type<scalar_coordinate_type>;;
    using real_coordinate_type = vector_type<float_type>;
    struct hash_functor
    {
        std::size_t operator()(Key const& n) const noexcept
        {
            return std::hash<value_type>()(n.index_);
        }
    };

public:  //static members

    static value_type compute_index(const coordinate_type& _c)noexcept {
        value_type idx =split_bits(static_cast<value_type>(_c.x()));
        for(int d=1;d<Dim;++d) 
        {
            idx |=(split_bits(static_cast<value_type>(_c[d])) << d);
        }
        return idx;
    }

public: 
    static coordinate_type coordinate(const value_type& code ) noexcept{
        coordinate_type c(0);
        for(int d=0;d<Dim;++d) c[d]=compress_bits(code >>d);
        return c;
    }

public: // Ctors

    Key() noexcept
    : index_(0) {}

    Key(value_type idx) noexcept
    : index_(idx) {}

    Key(int x, int y, int z) noexcept
    :Key(coordinate_type({x,y,z})) { }

    Key(coordinate_type x) noexcept
    : index_(compute_index(x)) {  }

    Key(const Key&) = default;
    Key(Key&&) = default;
    Key& operator=(const Key&) & = default;
    Key& operator=(Key&&) & = default;


public: //Access 

    coordinate_type coordinate()const noexcept {return coordinate(index_);}
    const value_type& index()const noexcept{return index_;}
    value_type& index() noexcept{return index_;}
    const value_type& getIndex()const noexcept {return index_;}
    value_type& getIndex() noexcept{return index_;}

    Key neighbor(const coordinate_type& _offset) const noexcept
    {
        const auto c =coordinate()+_offset;
        return Key(c);
    }
    
    Key child(int i) const noexcept
    {
        return (index_<<Dim)| static_cast<value_type>(i);
    }

    Key& operator++() noexcept { ++index_;return *this; }
    Key& operator--() noexcept { ++index_;return *this; }
    Key& operator+=(int) noexcept { ++index_; return *this;}
    Key& operator-=(int) noexcept { ++index_; return *this;}
    friend std::ostream& operator<<(std::ostream& os, const Key& _k)
    {
        os<<std::bitset<64>(_k.index_)<<" = "<<_k.index_
                 <<" coord =( "<<_k.coordinate()<<" )";
        return os;
    }


private: //Static 

    template<int nDim =Dim>
    static value_type split_bits(value_type w)
    {
        using tag=std::integral_constant<int,3>*;
        return split_bits_impl(w, tag(0));
    }
    template<int nDim =Dim>
    static value_type compress_bits(value_type w)
    {
        using tag=std::integral_constant<int,3>*;
        return compress_bits_impl(w, tag(0));
    }

    static value_type 
    split_bits_impl(value_type w, std::integral_constant<int, 3>*)  
    noexcept
    {
        w &=                0x00000000001fffff; 
        w = (w | w << 32) & 0x001f00000000ffff;  
        w = (w | w << 16) & 0x001f0000ff0000ff;  
        w = (w | w <<  8) & 0x010f00f00f00f00f; 
        w = (w | w <<  4) & 0x10c30c30c30c30c3; 
        w = (w | w <<  2) & 0x1249249249249249;
        return w;
    }
    static scalar_coordinate_type 
    compress_bits_impl(value_type w,std::integral_constant<int, 3>*) noexcept
    {
        w &=                  0x1249249249249249;
        w = (w ^ (w >> 2))  & 0x30c30c30c30c30c3;
        w = (w ^ (w >> 4))  & 0xf00f00f00f00f00f;
        w = (w ^ (w >> 8))  & 0x00ff0000ff0000ff;
        w = (w ^ (w >> 16)) & 0x00ff00000000ffff;
        w = (w ^ (w >> 32)) & 0x00000000001fffff;
        return static_cast<scalar_coordinate_type>(w);
    }
    static value_type 
    split_bits_impl(value_type x, std::integral_constant<int, 2>*)  
    noexcept
    {
        x = (x | (x << 16)) & 0x0000FFFF0000FFFF;
        x = (x | (x << 8))  & 0x00FF00FF00FF00FF;   
        x = (x | (x << 4))  & 0x0F0F0F0F0F0F0F0F; 
        x = (x | (x << 2))  & 0x3333333333333333;
        x = (x | (x << 1))  & 0x5555555555555555;
        return x;
    }

    static scalar_coordinate_type 
    compress_bits_impl(value_type w,std::integral_constant<int, 2>*) noexcept
    {
		w &=                  0x5555555555555555;
		w = (w ^ (w >> 1))  & 0x3333333333333333;
		w = (w ^ (w >> 2))  & 0x0f0f0f0f0f0f0f0f;
		w = (w ^ (w >> 4))  & 0x00ff00ff00ff00ff;
		w = (w ^ (w >> 8))  & 0x0000ffff0000ffff;
		w = (w ^ (w >> 16)) & 0x00000000ffffffff;
        return static_cast<scalar_coordinate_type>(w);
    }

private:
    value_type index_;
};

// binary operators
template<int Dim>
Key<Dim> operator+(int n, Key<Dim> k) noexcept
{ return k+=n;}

template<int Dim>
Key<Dim> operator-(int n, Key<Dim> k) noexcept
{ return k-=n; }

// relational operators
template<int Dim>
constexpr bool operator==(const Key<Dim>& l, const Key<Dim>& r) noexcept
{ return l.index() == r.index(); }

template<int Dim>
constexpr bool operator!=(const Key<Dim>& l, const Key<Dim>& r) noexcept
{ return l.index() != r.index(); }

template<int Dim>
constexpr bool operator<(const Key<Dim>& l, const Key<Dim>& r) noexcept
{ return l.index() < r.index(); }

template<int Dim>
constexpr bool operator<=(const Key<Dim>& l, const Key<Dim>& r) noexcept
{ return l.index() <= r.index(); }

template<int Dim>
constexpr bool operator>(const Key<Dim>& l, const Key<Dim>& r) noexcept
{ return l.index() > r.index(); }

template<int Dim>
constexpr bool operator>=(const Key<Dim>& l, const Key<Dim>& r) noexcept
{ return l.index() >= r.index(); }

}

#endif

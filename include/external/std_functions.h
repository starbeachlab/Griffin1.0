/////////////////////////////////////////////////////////////////////////
//!									 
//  This file is part of the GRIFFIN program				 
//									 
//  Copyright (C) 2010 by Rene Staritzbichler		      	 
//  rene.staritzbichler@biophys.mpg.de			       	 
//									 
//  GRIFFIN is free software; you can redistribute it and/or modify	 
//  it under the terms of the GNU General Public License as published by 
//  the Free Software Foundation; either version 2 of the License, or 
//  (at your option) any later version.				 
//									 
//  GRIFFIN is distributed in the hope that it will be useful,	 
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	 
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	 
//  GNU General Public License for more details.			 
//!									 
//!									 
//!									 
//!									 
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#ifndef EXTERNAL_STD_FUNCTIONS_H
#define EXTERNAL_STD_FUNCTIONS_H

#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <map>

#include <string.h>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"


template< typename t_DATA>
inline
std::vector< t_DATA> &operator *= ( std::vector< t_DATA> &VEC, const float &FACTOR)
{
    DebugWrite( __PRETTY_FUNCTION__);
    for( typename std::vector< t_DATA>::iterator itr = VEC.begin(); itr != VEC.end(); ++itr)
    {
        *itr *= FACTOR;
    }
    return VEC;
}

template< typename t_DATA>
inline
std::vector< t_DATA> operator * ( const std::vector< t_DATA> &VEC, const float &FACTOR)
{
    DebugWrite( __PRETTY_FUNCTION__);
    std::vector< t_DATA> vec( VEC.size());
    typename std::vector< t_DATA>::const_iterator o_itr = VEC.begin();
    for( typename std::vector< t_DATA>::iterator itr = vec.begin(); itr != vec.end(); ++itr, ++o_itr)
    {
        *itr = FACTOR * *o_itr;
    }
    return vec;
}

template< typename t_DATA>
inline
std::vector< t_DATA> operator * ( const float &FACTOR, const std::vector< t_DATA> &VEC)
{
    return VEC * FACTOR;
}


template< typename t_DATA>
inline
std::string GetStdVectorClassName( const std::vector< t_DATA> &VEC)
{ return "std::vector";}


template< typename t_DATA>
inline
std::istream& operator >> ( std::istream &STREAM, std::vector< t_DATA> &VEC)
{
  std::string str;
  STREAM >> str;
  assert( str == GetStdVectorClassName< t_DATA>( VEC));
  size_t size;
  STREAM >> size;
  VEC.clear();
  VEC.resize( size);
  for( typename std::vector< t_DATA>::iterator itr = VEC.begin(); itr != VEC.end(); ++itr)
  {
	  STREAM >> *itr;
  }
  return STREAM;
}


template< typename t_DATA>
inline
std::ostream& operator << ( std::ostream &STREAM, const std::vector< t_DATA> &VEC)
{
  STREAM << "std::vector" << std::endl;
  STREAM << VEC.size() << std::endl;
  for( typename std::vector< t_DATA>::const_iterator itr( VEC.begin()); itr != VEC.end(); ++itr)
  {
	  STREAM << *itr << "   ";
  }
  STREAM << std::endl;
  return STREAM;
}


template< typename t_FIRST, typename t_SECOND>
inline
std::ostream& operator << ( std::ostream &STREAM, const std::pair< t_FIRST, t_SECOND> &PAIR)
{
    STREAM << "std::pair" << std::endl;
    STREAM << PAIR.first << "   " << PAIR.second << std::endl;
    return STREAM;
}

template< typename t_FIRST, typename t_SECOND>
inline
std::istream& operator >> ( std::istream &STREAM, std::pair< t_FIRST, t_SECOND> &PAIR)
{
	std::string str;
	STREAM >> str;
	if( str != "std::pair")
    {
		std::cout << "====> expected \'std::pair\' as id, but got <" << str << "> " << std::endl;
		exit( -1);
    }
    STREAM >> PAIR.first >> PAIR.second;
    return STREAM;
}

template< typename t_KEY, typename t_VALUE>
inline
std::ostream& operator << ( std::ostream &STREAM, const std::map< t_KEY, t_VALUE> &MAP)
{
    STREAM << "std::map" << std::endl;
    STREAM << MAP.size() << std::endl;
    for( typename std::map< t_KEY, t_VALUE>::const_iterator itr = MAP.begin(); itr != MAP.end(); ++itr)
    {
        STREAM << itr->first << "   " << itr->second << std::endl;
    }
    return STREAM;
}

template< typename t_KEY, typename t_VALUE>
inline
std::istream& operator >> ( std::istream &STREAM, std::map< t_KEY, t_VALUE> &MAP)
{
	std::string
		str;
    int
		size;
    t_KEY
		key;
    t_VALUE
		value;

    STREAM >> str;

    if( str != "std::map")
    {
    	std::cout << "====> expected std::map, but got: " << str << std::endl;
    }

    STREAM >> size;

    for( int i = 0; i < size; ++i)
    {
    	STREAM >> key >> value;
    	MAP.insert( std::make_pair( key, value));
    }
    return STREAM;
}

template< typename t_KEY, typename t_VALUE>
inline
std::ostream& operator << ( std::ostream &STREAM, const std::multimap< t_KEY, t_VALUE> &MAP)
{
    STREAM << "std::multimap" << std::endl;
    for( typename std::multimap< t_KEY, t_VALUE>::const_iterator itr = MAP.begin(); itr != MAP.end(); ++itr)
    {
        STREAM << itr->first << "   " << itr->second << std::endl;
    }
    return STREAM;
}

inline
std::istream& operator >> ( std::istream &STREAM, std::vector< float> &VEC)
{
  std::string str;
  STREAM >> str;
  if( str == "std::vector<")
  {
      STREAM >> str;
      assert( str == "float>"); 
  }
  else if( str != "std::vector")
  {
      std::cerr << __FUNCTION__ << " expected std::vector, got: <" << str << ">" << std::endl;
      exit( -1);
  }
  size_t size;
  STREAM >> size;
  std::vector< float> tmp( size);
  for( size_t i = 0; i < size; ++i)
    { STREAM >> tmp[i];}
  VEC = tmp;
  return STREAM;
}

inline
std::ostream& operator << ( std::ostream &STREAM, const std::vector< float> &VEC)
{
  STREAM << "std::vector" << std::endl;
  STREAM << VEC.size() << std::endl;
  for( std::vector< float>::const_iterator itr = VEC.begin(); itr != VEC.end(); ++itr)
  {
	  STREAM << *itr << "   ";
  }
  STREAM << std::endl;
  return STREAM;
}


inline
std::istream& operator >> ( std::istream &STREAM, std::vector< size_t> &VEC)
{
  std::string str;
  STREAM >> str;
  if( str == "std::vector<")
  {
      STREAM >> str;
      assert( str == "size_t>");
  }
  else if( str != "std::vector")
  {
      std::cerr << __FUNCTION__ << " expected std::vector, got <" << str << ">" << std::endl;
      exit( -1);
  }
  size_t size;
  STREAM >> size;
  std::vector< size_t> tmp( size);
  for( size_t i = 0; i < size; ++i)
    { STREAM >> tmp[i];}
  VEC = tmp;
  return STREAM;
}

inline
std::ostream& operator << ( std::ostream &STREAM, const std::vector< size_t> &VEC)
{
  STREAM << "std::vector" << std::endl;
  STREAM << VEC.size() << std::endl;
  for( std::vector< size_t>::const_iterator itr = VEC.begin(); itr != VEC.end(); ++itr)
  {
	  STREAM << *itr << "   ";
  }
  STREAM << std::endl;
  return STREAM;
}


template< typename t_TYPE>
inline
std::vector< t_TYPE> Glue( const std::vector< t_TYPE> &V1, const std::vector< t_TYPE> &V2)
{
    DebugWrite( __PRETTY_FUNCTION__);
    std::vector< t_TYPE> vec( V1);
    vec.insert( vec.end(), V2.begin(), V2.end());
    return vec;
}


template< typename t_FIRST, typename t_SECOND>
inline
std::pair< t_FIRST, t_SECOND>  operator += ( const std::pair< t_FIRST, t_SECOND> &PAIR_A, const std::pair< t_FIRST, t_SECOND> &PAIR_B)
{
    DebugWrite( __PRETTY_FUNCTION__);
    std::pair< t_FIRST, t_SECOND> pair;
    pair.first = PAIR_A.first + PAIR_B.first;
    pair.second = PAIR_B.second + PAIR_B.second;
    return pair;
}

template< typename t_FIRST, typename t_SECOND>
inline
std::vector< t_SECOND> CastVector( const std::vector< t_FIRST> &SECOND)
{
    DebugWrite( __PRETTY_FUNCTION__);
    std::vector< t_SECOND> vector( SECOND.size());
    typename std::vector< t_SECOND>::iterator first_itr( vector.begin());
    typename std::vector< t_FIRST>::const_iterator second_itr( SECOND.begin());
    for( ; second_itr != SECOND.end(); ++first_itr, ++second_itr)
    {
        *first_itr = t_SECOND( **second_itr);
    }
    return vector;
}

template< typename t_FIRST, typename t_SECOND>
inline
std::map< t_FIRST, t_SECOND> &operator /= ( std::map< t_FIRST, t_SECOND> &MAP, const t_SECOND &OBJECT)
{
    DebugWrite( __PRETTY_FUNCTION__);
    for( typename std::map< t_FIRST, t_SECOND>::iterator itr = MAP.begin(); itr != MAP.end(); ++itr)
    {
        itr->second /= OBJECT;
    }
    return MAP;
}

template< typename t_TYPE>
inline
std::vector< t_TYPE> Merge( const std::vector< t_TYPE> &V1, const std::vector< t_TYPE> &V2)
{
    DebugWrite( __PRETTY_FUNCTION__);
    std::vector< t_TYPE> v( V1);
    v.insert( v.end(), V2.begin(), V2.end());
    return v;
}

template< typename t_TYPE>
inline
std::vector< t_TYPE> RemoveDuplicates( const std::vector< t_TYPE> &VECTOR)
{
    DebugWrite( __PRETTY_FUNCTION__);
    std::vector< t_TYPE> v( VECTOR), result;
    typename std::vector< t_TYPE>::iterator new_end = unique(v.begin(), v.end());
    result.insert( result.end(), v.begin(), new_end);
    return result;
}

template< typename t_TYPE>
inline
bool IsElement( const std::vector< t_TYPE> &VECTOR, const t_TYPE &VALUE)
{
    return(  std::find( VECTOR.begin(), VECTOR.end(), VALUE) != VECTOR.end());
}

template< typename t_FIRST, typename t_SECOND>
struct LargerThanForPairs
{
    bool operator()( const std::pair< t_FIRST, t_SECOND> &FIRST, const std::pair< t_FIRST, t_SECOND> &SECOND)
    {
        if( FIRST.first > SECOND.first)
        {
            return true;
        }
        else if( FIRST.first == SECOND.first && FIRST.second > SECOND.second)
        {
            return true;
        }
        return false;
    }
};

template< typename t_FIRST, typename t_SECOND>
struct SmallerThanForPairs
{
    bool operator()( const std::pair< t_FIRST, t_SECOND> &FIRST, const std::pair< t_FIRST, t_SECOND> &SECOND)
    {
        if( FIRST.first < SECOND.first)
        {
            return true;
        }
        else if( FIRST.first == SECOND.first && FIRST.second < SECOND.second)
        {
            return true;
        }
        return false;
    }
};


struct LargerThanForCharArrays
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};



#endif

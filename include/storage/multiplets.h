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


#ifndef MULTIPLETS_H_
#define MULTIPLETS_H_


template< typename T1, typename T2, typename T3>
struct Triplet
{
    T1 first;
    T2 second;
    T3 third;

    Triplet( const T1 &X, const T2 &Y, const T3 &Z)
    : first( X), second( Y), third( Z)
      {}
};

template< typename T1, typename T2, typename T3, typename T4>
struct Quartet
{
    T1 first;
    T2 second;
    T3 third;
    T4 forth;
};

template< typename T1, typename T2, typename T3>
inline
std::ostream &operator << ( std::ostream &STREAM, const Triplet< T1, T2, T3> &TRIPLET)
{
    STREAM << "(" << TRIPLET.first << ", " << TRIPLET.second << ", " << TRIPLET.third << ")";
    return STREAM;
}

template< typename T1, typename T2, typename T3, typename T4>
inline
std::ostream &operator << ( std::ostream &STREAM, const Quartet< T1, T2, T3, T4> &QUARTET)
{
    STREAM << "(" << QUARTET.first << ", " << QUARTET.second << ", " << QUARTET.third << ", " << QUARTET.forth << ")";
    return STREAM;
}



#endif /* MULTIPLETS_H_ */

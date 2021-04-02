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


#ifndef MATH_DOUBLE_FUNCTIONS_H
#define MATH_DOUBLE_FUNCTIONS_H


#include <vector>

#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "constants.h"


namespace math
{


  inline float Unsign( const float &X)
  {
    if( X < 0){ return -X;}
    return X;
  }

  //! rounds positive (!) float values, returns an size_t
  template< typename t_TYPE>
  inline t_TYPE Round( const float &X)
  {
    if( X < 0.0)
    { return t_TYPE( int( X - 0.5));}
    return t_TYPE( X + 0.5);
  }

  //! rounds float values, returns a float
  inline float Round( const float &X)
  {
    if( X < 0.0)
    { return float( int( X - 0.5));}
    return float( int( X + 0.5));
  }


  //! truncates a float value, returns an float
  inline float Truncate( const float &X)
  { return float( int( X));}

  //! random [MIN,MAX)
  inline float BasicRandom( const float &MIN, const float &MAX)
  {
      return float( rand() / ( float( RAND_MAX) + 1.0) * ( MAX - MIN) + MIN);
  }


  inline
  float Square( const float &VALUE)
    { return VALUE * VALUE;}


  template< typename T>
  std::vector< T> Devisors( const T &VALUE)
    { /* me no smart */
      std::vector< T> devisors;
      for( T i( 1); i <= 0.5 * VALUE; ++i) // improve limits!
    if( VALUE % i == 0)
      { devisors.push_back( i);}
      return devisors;
    }

  inline
  float AngleFromLawOfCosinus( const float &OPPOSITE_SIDE_LENGTH, const float &FIRST_NEXT_SIDE_LENGTH, const float &SECOND_NEXT_SIDE_LENGTH)
  { return acos( ( FIRST_NEXT_SIDE_LENGTH * FIRST_NEXT_SIDE_LENGTH + SECOND_NEXT_SIDE_LENGTH * SECOND_NEXT_SIDE_LENGTH - OPPOSITE_SIDE_LENGTH * OPPOSITE_SIDE_LENGTH) / ( 2.0 * FIRST_NEXT_SIDE_LENGTH * SECOND_NEXT_SIDE_LENGTH));}

  inline
  float CosinusAngleFromLawOfCosinus( const float &OPPOSITE_SIDE_LENGTH, const float &FIRST_NEXT_SIDE_LENGTH, const float &SECOND_NEXT_SIDE_LENGTH)
  { return ( FIRST_NEXT_SIDE_LENGTH * FIRST_NEXT_SIDE_LENGTH + SECOND_NEXT_SIDE_LENGTH * SECOND_NEXT_SIDE_LENGTH - OPPOSITE_SIDE_LENGTH * OPPOSITE_SIDE_LENGTH) / ( 2.0 * FIRST_NEXT_SIDE_LENGTH * SECOND_NEXT_SIDE_LENGTH);}

  inline
  float LawOfCosinus( const float &AN_1, const float &AN_2, const float &ANGLE)
  { return ( AN_1 * AN_1 + AN_2 * AN_2 - 2.0 * AN_1 * AN_2 * cos( ANGLE));}

  inline
  float ProjectedLawOfCosinus( const float &FIRST_SIDE_LENGTH, const float &FIRST_ANGLE, const float &SECOND_SIDE_LENGTH, const float &SECOND_ANGLE)
  { return FIRST_SIDE_LENGTH * cos( SECOND_ANGLE) + SECOND_SIDE_LENGTH * cos( FIRST_ANGLE);}

  inline
  float LengthFromLawOfSinus( const float &OPPOSITE_ANGLE, const float &OTHER_SIDE_LENGTH, const float &OTHER_OPPOSITE_ANGLE)
  { return OTHER_SIDE_LENGTH * sin( OPPOSITE_ANGLE) / sin( OTHER_OPPOSITE_ANGLE);}

  inline
  float AngleFromLawOfSinus( const float & OPPOSITE_SIDE_LENGTH, const float &OTHER_SIDE_LENGTH, const float &OTHER_ANGLE)
  { return asin( OPPOSITE_SIDE_LENGTH / OTHER_SIDE_LENGTH * sin( OTHER_ANGLE));}

  inline
  float RadiansToDegrees( const float &RAD)
  { return RAD * 180.0 / math::Pi;}

  inline
  float DegreesToRadians( const float &DEG)
  { return DEG * math::Pi / 180.0;}

  inline
  bool IsLarger( const float& X, const float &Y)
  { return X > Y;}

  inline
  bool IsSmaller( const float& X, const float &Y)
  { return X < Y;}

  inline
  bool IsEqual( const float& X, const float &Y)
  { return X == Y;}

  inline
  bool IsEqualWithinThreshold( const float& X, const float &Y, const float &THRESHOLD)
  { return fabs( X - Y) < THRESHOLD;}

  inline
  bool IsNonEqual( const float& X, const float &Y)
  { return X != Y;}

  inline
  float
  Modulo( const float &VALUE, const float &NOMINATOR)
  {
    return VALUE - float( int( VALUE / NOMINATOR)) * NOMINATOR;
  }

  inline
  float
  GaussDistribution( const float &X, const float &DELTA)
  {
	  return  1.0 / ( DELTA * sqrt( 2.0 * math::Pi)) * exp( -0.5 * ( X * X / ( DELTA * DELTA)));
  }

} // end namespace math


#endif

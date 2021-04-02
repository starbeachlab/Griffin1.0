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


#ifndef VECTOR_FUNCTIONS_H
#define VECTOR_FUNCTIONS_H


#include "../../include/math/vector.h"

namespace math
{
  //! returns normalized vector, conserves argument vector
  Vector NormalizedVector( const Vector &V);


  //! returns the scalar product of two vectors V1 and V2
  float ScalarProduct( const Vector &V1, const Vector &V2);

  //! write vector into std::ostream
  std::ostream &operator << ( std::ostream & STREAM, const Vector& V);

  //! read vector from std::istream
  std::istream &operator >> ( std::istream &STREAM, Vector &V);

  //! difference of two vectors
  Vector operator - ( const Vector &V1, const Vector &V2);

  //! sum of two vectors
  Vector operator + ( const Vector &V1, const Vector &V2);

  //! scalar product of two vectors
  float operator * ( const Vector &V1, const Vector &V2);

  //! multiply vector with factor
  Vector operator * ( const float &VALUE, const Vector &V);

  //! multiply vector with factor
  Vector operator * ( const Vector &V, const float &VALUE);

  //! devide vector by factor
  Vector operator / ( const Vector &V, const float &VALUE);

  //! devide vector by vector
  Vector operator / ( const Vector &V1, const Vector &V2);

  //! the squared distance of two points
  float SquaredDistance( const Vector &POS_A, const Vector &POS_B);

  //! the distance of two points
   float Distance( const Vector &POS_A, const Vector &POS_B);

   //! the angle between two vectors
   float Angle( const Vector &DIRECTION_A, const Vector &DIRECTION_B);

   //! boolean function checking whether the distance of two points is smaller than threshold
   bool IsDistanceSmallerThan( const Vector &POSITION_A, const Vector &POSITION_B, const float &THRESHOLD);

   //! the length of the projection of vector A on vector B
   float ProjectionOnVector( const math::Vector &A, const math::Vector &B);

   //! returns elementwise minimum between two vectors A and B
   Vector Min( const Vector &A, const Vector &B);

    //! returns elementwise maximum between two vectors A and B
    Vector Max( const Vector &A, const Vector &B);


} // end namespace math


#endif

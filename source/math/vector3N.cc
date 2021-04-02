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


#include <cassert>
#include "../../include/math/vector3N.h"


namespace math
{

  Vector3N::Vector3N()
    : Vector(3, 0.0)
  {}

  Vector3N::Vector3N( const Vector3N &V)
    : Vector( V)
  { assert( V.size() == 3);}

  Vector3N::Vector3N( const Vector &V)
    : Vector( V)
  { assert( V.size() == 3);}

  Vector3N::Vector3N( const std::vector< float> &V)
    : Vector( V)
  { assert( V.size() == 3);}

  Vector3N::Vector3N( const float &X)
    : Vector( 3)
  {
    ( *this)[0] = X;
    ( *this)[1] = X;
    ( *this)[2] = X;
  }

  Vector3N::Vector3N( const float &X, const float &Y, const float &Z)
    : Vector( 3)
  {
    ( *this)[0] = X;
    ( *this)[1] = Y;
    ( *this)[2] = Z;
  }

  Vector3N::~Vector3N(){}

  Vector3N *Vector3N::Clone() const
  { return new Vector3N( *this);}

  float &Vector3N::operator()( const size_t &ID)
  {
//    assert( ID < 3);
    return ( *this)[ ID];
  }

  const float& Vector3N::operator() ( const size_t &ID) const
  {
//    assert( ID < 3);
    return ( *this)[ ID];
  }

  Vector3N &Vector3N::operator = ( const std::vector< float> &VEC)
  {
      assert( VEC.size() == 3);
        ( *this)[0] = VEC[0];
        ( *this)[1] = VEC[1];
        ( *this)[2] = VEC[2];
        return *this;
  }


  std::istream &Vector3N::Read( std::istream &STREAM)
  { return Vector::Read( STREAM);}


  std::ostream &Vector3N::Write( std::ostream &STREAM) const
  { return Vector::Write( STREAM);}

} // end namespace math

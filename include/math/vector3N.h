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


#ifndef MATH_VECTOR3N_H
#define MATH_VECTOR3N_H

#include "vector.h"
#include "../../include/math/vector_functions.h"

namespace math
{

  class Vector3N
    : public Vector
  {
  public:
    Vector3N();

    Vector3N( const Vector3N &V);

    Vector3N( const Vector &V);
    Vector3N( const std::vector< float> &V);

    Vector3N( const float &X);

    Vector3N( const float &X, const float &Y, const float &Z);

    virtual ~Vector3N();

    virtual Vector3N *Clone() const;

    virtual const float& operator() ( const size_t &ID) const;

    virtual float& operator() ( const size_t &ID);

    Vector3N &operator = ( const std::vector< float> &VEC);

    virtual std::ostream &Write( std::ostream &STREAM) const;

    virtual std::istream &Read( std::istream &STREAM);

    //! insert all elements of Vector in this
    virtual Vector3N &Merge( const Vector3N &V)
    {
        std::cout << "calling this function of a vector of fixed size does not make any sense!" << std::endl;
        exit( -1);
        return *this;
    }

    //! remove multiple copies
    virtual Vector3N &Unique()
    {
        std::cout << "calling this function of a vector of fixed size does not make any sense!" << std::endl;
        exit( -1);
        return *this;
    }



  }; // class Vector3N

    inline
    Vector3N CrossProduct( const Vector3N &X, const Vector3N &Y)
    {
        return Vector3N
        (
                X( 1) * Y( 2) - X( 2) * Y( 1),
             X( 2) * Y( 0) - X( 0) * Y( 2),
             X( 0) * Y( 1) - X( 1) * Y( 0)
        );
    }

    //! Spatprodukt - the volume created from the three vectors
    inline
    float TripleProduct( const Vector3N &X, const Vector3N &Y, const Vector3N &Z)
    {
        return ScalarProduct( CrossProduct( X, Y), Z);
    }


} // namespace math

#endif

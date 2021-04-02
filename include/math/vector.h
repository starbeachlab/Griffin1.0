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


#ifndef MATH_VECTOR_H
#define MATH_VECTOR_H

#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include <boost/shared_ptr.hpp>

// forward declaration
namespace math
{ class Vector;}

//#include "../../include/string/io_string_functions.h"
//#include "../../include/math/vector_functions.h"
#include "../readwrite/stream_operator.h"
#include "../readwrite/class_name.h"
#include "double_functions.h"

namespace math
{

  ///////////////////////////////////////////////
  //! Wrapper class around std::vector< float>
  //! with secure element access.
  //!
  //! @author Rene Staritzbichler
  //! @date 19.11.2008
  //! @example "../../example/math/vector.cpp", "../../example/math/vector3N.cpp","../../example/math/vector_functions.cpp"
  ///////////////////////////////////////////////

  class Vector
    : public std::vector< float>, public StreamOperator, public readwrite::ClassName
  {
  public:
    //! default constructor
    Vector()
      : std::vector<float>(),
      StreamOperator(),
      ClassName()
      {}

    //! construction from size
    Vector( const size_t &SIZE)
  : std::vector< float>( SIZE),
    StreamOperator(),
    ClassName()
  { /*SetAll( 0.0);*/}

    //! construction from size and value
    Vector( const size_t &SIZE, const float &VALUE)
  : std::vector< float>( SIZE, VALUE),
    StreamOperator(),
    ClassName()
  { /*SetAll( 0.0);*/}

      Vector( const std::vector< size_t> &V)
        : std::vector< float>( V.size()),
        StreamOperator(),
        ClassName()
        {
          std::vector< size_t>::const_iterator v_itr( V.begin());
          for( std::vector< float>::iterator this_itr( this->begin()); this_itr != this->end() && v_itr != V.end(); ++v_itr, ++this_itr)
        { *this_itr = float( *v_itr);}
        }

      Vector( const std::vector< float> &V)
        : std::vector< float>( V),
        StreamOperator(),
        ClassName()
        {}


      //! copy constructor
        Vector( const Vector &V)
          : std::vector< float>( V),
          StreamOperator(),
          ClassName()
          {}


      //! virtual destructor
      virtual ~Vector(){}

      //! mutable element access
      virtual float& operator()( const size_t &ID);

      //! unmutable element access
      virtual const float& operator()( const size_t &ID) const;

      //! add another vector
      virtual Vector& operator += ( const Vector &V);

      //! substract another vector
      virtual Vector& operator -= ( const Vector &V);

//      //! set this vector to the values of another vector
//      virtual Vector operator = ( const Vector &V);

      //! multiply each element of this vector with VALUE
      virtual Vector& operator *= ( const float &VALUE);
      //! multiply each element of this vector with VALUE
 //     virtual Vector operator *= ( const float &VALUE) const;
      //! devide each element of this vector by VALUE
      virtual Vector& operator /= ( const float &VALUE);
      //! substract VALUE from each element of this vector
      virtual Vector& operator -= ( const float &VALUE);
      //! add VALUE to each element of this vector
      virtual Vector& operator += ( const float &VALUE);
      //! set each element of this vector to VALUE
      virtual Vector &operator = ( const float &VALUE);
//      //! check whether this vector is equal to V
//      virtual bool operator == ( const Vector &V);

      //! checks whether condition '<' is fullilled for all elements
      virtual bool operator < ( const Vector &V) const;

      //! checks whether condition '<=' is fullilled for all elements
      virtual bool operator <= ( const Vector &V) const;

      //! checks whether condition ">' is fullilled for all elements
      virtual bool operator > ( const Vector &V) const;

      //! checks whether condition '>=' is fullilled for all elements
      virtual bool operator >= ( const Vector &V) const;

      virtual float SumOfElements() const;

      virtual float SubSumOfElements( const size_t &BEGIN, const size_t &NR_ELEMENTS) const;

      virtual float ProductOfElements() const;

      virtual float SubProductOfElements( const size_t &BEGIN, const size_t &NR_ELEMENTS) const;

      virtual Vector SubVector( const size_t &BEGIN, const size_t &NR_ELEMENTS) const;

      virtual float MeanOfElements() const;

      virtual std::pair< float, float> RmsdOfElements() const;

      virtual float MaxElement() const;

      virtual float MinElement() const;

      virtual float DeltaOfElements() const;

      virtual void SetAll( const float &VALUE);

      //! check whether length of vector is larger zero (faster than explicit calculation of the length)
      virtual bool IsLengthLargerZero() const;

      //! returns the length of the vector
      virtual float Length() const;

      //! returns the squared length of the vector
      virtual float SquaredLength() const;

      //! normalize the vector
      virtual Vector &Normalize();

      //! return normalized copy of vector
      virtual Vector NormalizedCopy() const;

      //! returns length and normalizes vector
      virtual float GetLengthAndNormalize();

      //! write vector to ostream
      virtual std::ostream& Write( std::ostream &STREAM) const;

      //! read vector from istream
      virtual std::istream& Read( std::istream &STREAM);

      //! randomize all elements within limits
      virtual Vector &Randomize( const float &MIN, const float &MAX);

      //! randomize all elements within limits
      virtual Vector &Randomize( const Vector &MIN, const Vector &MAX);

      //! set vector length to MAX if larger MAX
      virtual Vector &SetToLength( const float &MAX);

      //! sets all negative elements to -1, all positive to +1, 0 remains 0
      virtual Vector &SetNonZeroElementsToMinusOrPlusOne();

      //! devide each element of this vector by element of VECTOR
      virtual Vector &operator /= ( const Vector &VECTOR);

      //! insert all elements of Vector in this
      virtual Vector &Merge( const Vector &V);

      //! checks whether all element in the vector are equal VALUE
      virtual bool AreAllElements( const float &VALUE) const;


//      //! remove multiple copies
//      virtual Vector &RemoveDuplicates();



  }; // end class Vector

} // end namespace math

#endif

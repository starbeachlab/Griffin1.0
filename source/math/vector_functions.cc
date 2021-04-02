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
#include "../../include/math/vector_functions.h"

namespace math
{
  Vector NormalizedVector( const Vector &V)
  {
    Vector v( V);
    return v.Normalize();
  }

  float ScalarProduct( const Vector &V1, const Vector &V2)
  {
    assert( V1.size() == V2.size());
    float scalar(0.0);
    for( std::vector< float>::const_iterator itr1( V1.begin()), itr2( V2.begin()); itr1 != V1.end() && itr2 != V2.end(); ++itr1, ++itr2)
      { scalar += ( *itr1) * ( *itr2);}
    return scalar;
  }

  std::ostream &operator << ( std::ostream & STREAM, const Vector& V)
  { return V.Write( STREAM);}

  std::istream &operator >> ( std::istream &STREAM, Vector &V)
  { return V.Read( STREAM);}

  Vector operator - ( const Vector &V1, const Vector &V2)
    {
      assert( V1.size() == V2.size());
      Vector difference( V1.size());
      std::vector< float>::iterator itr_diff( difference.begin());
      for( std::vector<float>::const_iterator itr_v1( V1.begin()), itr_v2( V2.begin()); itr_v1 != V1.end(); ++itr_v1, ++itr_v2, ++itr_diff)
      { *itr_diff = *itr_v1 - *itr_v2;}
      return difference;
  }

  Vector operator + ( const Vector &V1, const Vector &V2)
    {
      assert( V1.size() == V2.size());
      Vector sum( V1.size());
      std::vector< float>::iterator itr_sum( sum.begin());
      for( std::vector<float>::const_iterator itr_v1( V1.begin()), itr_v2( V2.begin()); itr_v1 != V1.end(); ++itr_v1, ++itr_v2, ++itr_sum)
      { *itr_sum = *itr_v1 + *itr_v2;}
      return sum;
  }

  float operator * ( const Vector &V1, const Vector &V2)
    {
      assert( V1.size() == V2.size());
      float scalar( 0.0);
      for( std::vector< float>::const_iterator itr_v1( V1.begin()), itr_v2( V2.begin()); itr_v1 != V1.end(); ++itr_v1, ++itr_v2)
      { scalar += *itr_v1 * *itr_v2;}
      return scalar;
    }

  //! multiply vector with factor
  Vector operator * ( const float &VALUE, const Vector &V)
  {
      Vector vec( V);
      for( std::vector< float>::iterator itr = vec.begin(); itr != vec.end(); ++itr)
      {
          *itr *= VALUE;
      }
      return vec;
  }

  //! multiply vector with factor
  Vector operator * ( const Vector &V, const float &VALUE)
{
      Vector vec( V);
      for( std::vector< float>::iterator itr = vec.begin(); itr != vec.end(); ++itr)
      {
          *itr *= VALUE;
      }
      return vec;
}

  Vector operator / ( const Vector &V, const float &VALUE)
  {
	  Vector vec( V);
      for( std::vector< float>::iterator itr = vec.begin(); itr != vec.end(); ++itr)
      {
          *itr /= VALUE;
      }
      return vec;
  }


  Vector operator / ( const Vector &V1, const Vector &V2)
  {
      assert( V1.size() == V2.size());
      Vector vec( V1.size());
      std::vector< float>::iterator itr = vec.begin();
      for( std::vector< float>::const_iterator itr_v1( V1.begin()), itr_v2( V2.begin()); itr_v1 != V1.end(); ++itr_v1, ++itr_v2, ++itr)
      { *itr = *itr_v1 / *itr_v2;}
      return vec;
  }

  float SquaredDistance( const Vector &POSITION_A, const Vector &POSITION_B)
  {
    assert( POSITION_A.size() == POSITION_B.size());
    float dist( 0.0);
    for( std::vector< float>::const_iterator itr_a( POSITION_A.begin()), itr_b( POSITION_B.begin()); itr_a != POSITION_A.end() && itr_b != POSITION_B.end(); ++itr_a, ++itr_b)
      { dist += ( *itr_b - *itr_a) * ( *itr_b - *itr_a);}
    return dist;
  }

  float Distance( const Vector &POSITION_A, const Vector &POSITION_B)
  {
    assert( POSITION_A.size() == POSITION_B.size());
    float dist( 0.0);
    for( std::vector< float>::const_iterator itr_a( POSITION_A.begin()), itr_b( POSITION_B.begin()); itr_a != POSITION_A.end() && itr_b != POSITION_B.end(); ++itr_a, ++itr_b)
      { dist += ( *itr_b - *itr_a) * ( *itr_b - *itr_a);}
    return sqrt( dist);
  }

  float Angle( const Vector &DIRECTION_A, const Vector &DIRECTION_B)
  {
      float 
	  a_sqr =  DIRECTION_A.SquaredLength(),
	  b_sqr =  DIRECTION_B.SquaredLength();

      if( a_sqr == 0 || b_sqr == 0)
      {
          std::cout << "====> angle of vectors of size zero requested!" << std::endl;
          return std::numeric_limits< float>::max(); // TODO: return undefined
      }
      float scalar( ScalarProduct( DIRECTION_A, DIRECTION_B) / sqrt( a_sqr * b_sqr));
      if( scalar > 1.0 && scalar < 1.00001)
      {
    	  scalar = 1.0;
      }
      else if( scalar < -1.0 && scalar > -1.00001)
      {
          scalar = -1.0;
      }
      else if( scalar > 1.0 || scalar < -1.0)
      {
          std::cout << "===> undefined values in scalar product: " << scalar << std::endl;
          exit( -1);
      }
    return acos( scalar);
  }

  bool IsDistanceSmallerThan( const Vector &POSITION_A, const Vector &POSITION_B, const float &THRESHOLD)
  {
        assert( POSITION_A.size() == POSITION_B.size());
        float squared_distance( 0.0);
        for( std::vector< float>::const_iterator itr_a( POSITION_A.begin()), itr_b( POSITION_B.begin()); itr_a != POSITION_A.end() && itr_b != POSITION_B.end(); ++itr_a, ++itr_b)
          { squared_distance += ( *itr_b - *itr_a) * ( *itr_b - *itr_a);}
        return squared_distance < THRESHOLD * THRESHOLD;
  }

  //! the length of the projection of vector A on vector B
  float ProjectionOnVector( const math::Vector &A, const math::Vector &B)
  {
       return ScalarProduct( A, B) / B.Length();
  }

  //! returns elementwise minimum between two vectors A and B
   Vector Min( const Vector &A, const Vector &B)
   {
	      assert( A.size() == B.size());
	      Vector vec( A.size());
	      std::vector< float>::iterator itr = vec.begin();
	      for( std::vector< float>::const_iterator itr_v1( A.begin()), itr_v2( B.begin()); itr_v1 != A.end(); ++itr_v1, ++itr_v2, ++itr)
	      { *itr = std::min( *itr_v1, *itr_v2);}
	      return vec;

   }

    //! returns elementwise maximum between two vectors A and B
    Vector Max( const Vector &A, const Vector &B)
    {
 	      assert( A.size() == B.size());
 	      Vector vec( A.size());
 	      std::vector< float>::iterator itr = vec.begin();
 	      for( std::vector< float>::const_iterator itr_v1( A.begin()), itr_v2( B.begin()); itr_v1 != A.end(); ++itr_v1, ++itr_v2, ++itr)
 	      { *itr = std::max( *itr_v1, *itr_v2);}
 	      return vec;

    }


} // end namespace math


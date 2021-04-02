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


#include "../../include/math/vector.h"


namespace math
{

  float& Vector::operator()( const size_t &ID)
  {
    assert( ID < this->size());
    return this->operator[]( ID);
  }


  const float& Vector::operator()( const size_t &ID) const
  {
    assert( ID < this->size());
    return this->operator[]( ID);
  }

  Vector &Vector::operator += ( const Vector &V)
  {
      assert( this->size() == V.size());
      std::vector< float>::const_iterator v_itr( V.begin());
      for( std::vector< float>::iterator itr = this->begin(); itr != this->end(); ++itr, ++v_itr)
      {
          *itr += *v_itr;
      }
      return *this;
  }


  Vector &Vector::operator -= ( const Vector &V)
  {
    assert( this->size() == V.size());
    std::vector< float>::const_iterator v_itr( V.begin());
    for( std::vector< float>::iterator itr = this->begin(); itr != this->end(); ++itr, ++v_itr)
    {
    	*itr -= *v_itr;
    }
    return *this;
  }


  Vector &Vector::operator *= ( const float &VALUE)
  {
      for( std::vector< float>::iterator itr = this->begin(); itr != this->end(); ++itr)
      {
          *itr *= VALUE;
      }
      return *this;
  }

//  //! multiply each element of this vector with VALUE
//  Vector Vector::operator *= ( const float &VALUE) const
//  {
//      Vector vec( *this);
//      for( std::vector< float>::iterator itr = vec.begin(); itr != vec.end(); ++itr)
//      {
//          *itr *= VALUE;
//      }
//      return vec;
//  }

  Vector &Vector::operator /= ( const float &VALUE)
  {
      float inverse( 1.0 / VALUE);
      for( std::vector< float>::iterator itr = this->begin(); itr != this->end(); ++itr)
      {
          *itr *= inverse;
      }
      return *this;
  }

  Vector &Vector::operator /= ( const Vector &VEC)
  {
      assert( this->size() == VEC.size());
      std::vector< float>::const_iterator vec_itr = VEC.begin();
      for( std::vector< float>::iterator this_itr = this->begin(); this_itr != this->end(); ++this_itr, ++vec_itr)
      {
          *this_itr /= *vec_itr;
      }
      return *this;
  }

  Vector &Vector::operator -= ( const float &VALUE)
  {
      for( std::vector< float>::iterator itr = this->begin(); itr != this->end(); ++itr)
      {
          *itr -= VALUE;
      }
      return *this;
  }

  Vector &Vector::operator += ( const float &VALUE)
  {
       for( std::vector< float>::iterator itr = this->begin(); itr != this->end(); ++itr)
       {
           *itr += VALUE;
       }
       return *this;
  }

  Vector &Vector::operator = ( const float &VALUE)
  {
    for( std::vector< float>::iterator itr = this->begin(); itr != this->end(); ++itr)
      {
    *itr = VALUE;
      }
    return *this;
  }

  bool Vector::operator < ( const Vector &V) const
  {
      assert( this->size() == V.size());
      for( std::vector< float>::const_iterator this_itr( this->begin()), v_itr( V.begin()); this_itr != this->end(); ++this_itr, ++v_itr)
          if( *this_itr >= *v_itr)
          { return false;}
      return true;
  }

  bool Vector::operator <= ( const Vector &V) const
  {
      assert( this->size() == V.size());
      for( std::vector< float>::const_iterator this_itr( this->begin()), v_itr( V.begin()); this_itr != this->end(); ++this_itr, ++v_itr)
          if( *this_itr > *v_itr)
          { return false;}
      return true;
  }

  bool Vector::operator > ( const Vector &V) const
  {
      assert( this->size() == V.size());
      for( std::vector< float>::const_iterator this_itr( this->begin()), v_itr( V.begin()); this_itr != this->end(); ++this_itr, ++v_itr)
          if( *this_itr <= *v_itr)
          { return false;}
      return true;
  }

  bool Vector::operator >= ( const Vector &V) const
  {
      assert( this->size() == V.size());
      for( std::vector< float>::const_iterator this_itr( this->begin()), v_itr( V.begin()); this_itr != this->end(); ++this_itr, ++v_itr)
          if( *this_itr < *v_itr)
          { return false;}
      return true;
  }


  float Vector::SumOfElements() const
  {
    float sum( 0.0);
    for( std::vector< float>::const_iterator itr( this->begin()); itr != this->end(); ++itr)
      { sum += *itr;}
    return sum;
  }

  float Vector::SubSumOfElements( const size_t &BEGIN, const size_t &NR_ELEMENTS) const
  {
    size_t count( 0);
    float sum( 0.0);
    for( std::vector< float>::const_iterator itr( this->begin() + BEGIN); itr != this->end() && count < NR_ELEMENTS; ++itr, ++count)
      { sum += *itr;}
    return sum;
  }


  float Vector::ProductOfElements() const
  {
    float product( 1.0);
    for( std::vector< float>::const_iterator itr( this->begin()); itr != this->end(); ++itr)
      { product *= *itr;}
    return product;
  }

  float Vector::SubProductOfElements( const size_t &BEGIN, const size_t &NR_ELEMENTS) const
  {
    size_t count( 0);
    float product( 1.0);
    for( std::vector< float>::const_iterator itr( this->begin() + BEGIN); itr != this->end() && count < NR_ELEMENTS; ++itr, ++count)
      { product *= *itr;}
    return product;
  }

  Vector Vector::SubVector( const size_t &BEGIN, const size_t &NR_ELEMENTS) const
  {
    Vector sub( NR_ELEMENTS);
    std::vector< float>::iterator new_itr( sub.begin());
    for( std::vector< float>::const_iterator this_itr( this->begin() + BEGIN); this_itr != this->end() && new_itr != sub.end(); ++this_itr, ++new_itr)
      { *new_itr = *this_itr;}
    return sub;
  }

  float Vector::MeanOfElements() const
  { return SumOfElements() / this->size();}


  std::pair< float, float> Vector::RmsdOfElements() const
  { // could be optimized
    float rmsd( 0.0);
    float mean( MeanOfElements());
    for( std::vector< float>::const_iterator itr = this->begin(); itr != this->end(); ++itr)
      { rmsd += math::Square( ( *itr) - mean);}
    return std::pair< float, float>( mean, sqrt( rmsd / float( this->size() - 1)));
  }

  float Vector::MaxElement() const
  {
    float max( *this->begin());
    for( std::vector< float>::const_iterator itr = this->begin() + 1; itr != this->end(); ++itr)
      if( *itr > max)
    { max = *itr;}
    return max;
  }

  float Vector::MinElement() const
  {
    float min( *this->begin());
    for( std::vector< float>::const_iterator itr = this->begin() + 1; itr != this->end(); ++itr)
      if( *itr < min)
    { min = *itr;}
    return min;
  }

  float Vector::DeltaOfElements() const
  { return MaxElement() - MinElement();}



  std::ostream& Vector::Write( std::ostream &STREAM) const
  {
	  STREAM << mystr::GetClassName( std::string( __PRETTY_FUNCTION__)) << std::endl;
	  STREAM << this->size() << std::endl;
	  for( std::vector< float>::const_iterator itr = this->begin(); itr != this->end(); ++itr)
      {
		  STREAM.width( 16);
		  STREAM << *itr << "   ";
      }
	  STREAM << std::endl;
	  return STREAM;
  }


  std::istream& Vector::Read( std::istream &STREAM)
  {
    std::string str;
    STREAM >> str;
    if( str != mystr::GetClassName( std::string( __PRETTY_FUNCTION__)))
    {
    	std::cout << "===> not the correct id <" << str << "> in " << __PRETTY_FUNCTION__ << std::endl;
    	STREAM >> str;
    	std::cout << "try next: " << str << std::endl;
    	exit( -1);
    }
    size_t size;
    STREAM >> size;
    *this = Vector( size);
    for( std::vector< float>::iterator itr = this->begin(); itr != this->end(); ++itr)
    {
    	STREAM >> *itr;
    }
    return STREAM;
  }


  bool Vector::IsLengthLargerZero() const
  {
	  for( std::vector< float>::const_iterator itr( this->begin()); itr != this->end(); ++itr)
		  if( *itr != 0)
		  { return true;}
	  return false;
  }


  float Vector::Length() const
  {
	  float sum = 0.0;
	  for( std::vector< float>::const_iterator itr = this->begin(); itr != this->end(); ++itr)
	  {
		  sum += *itr * *itr;
	  }
	  return sqrt( sum);
  }

  float Vector::SquaredLength() const
  {
	  float sum = 0.0;
	  for( std::vector< float>::const_iterator itr = this->begin(); itr != this->end(); ++itr)
	  {
		  sum += *itr * *itr;
	  }
	  return sum;
  }


  Vector &Vector::Normalize()
  {
    *this *= 1.0 / Length();
    return *this;
  }

  Vector Vector::NormalizedCopy() const
  {
      Vector copy( *this);
      copy *= 1.0 / Length();
      return copy;
  }

  float Vector::GetLengthAndNormalize()
  {
	  float length = Length();
      *this *= 1.0 / length;
      return length;
  }

  void Vector::SetAll( const float& X)
  {
      for( std::vector< float>::iterator itr = this->begin(); itr != this->end(); ++itr)
      { *itr = X;}
  }

  Vector &Vector::Randomize( const float &MIN, const float &MAX)
  {
      for( std::vector< float>::iterator itr = this->begin(); itr != this->end(); ++itr)
      {  *itr = math::BasicRandom( MIN, MAX);}
      return *this;
  }

  //! randomize all elements within limits
  Vector &Vector::Randomize( const Vector &MIN, const Vector &MAX)
  {
      std::vector< float>::const_iterator min_itr = MIN.begin(), max_itr = MAX.begin();
      for( std::vector< float>::iterator itr = this->begin(); itr != this->end(); ++itr, ++max_itr, ++min_itr)
      {
          *itr = math::BasicRandom( *min_itr, *max_itr);
      }
      return *this;
  }


  Vector &Vector::SetToLength( const float &MAX)
  {
      float length( Length());
      DebugWrite( "Vector::SetToLength(): length_before: " << length);
      if( length > 0 && MAX >= 0.0 && MAX < std::numeric_limits< float>::max())
      {
    	  *this *= ( MAX / length);
#ifdef DEBUG
    	  std::cout << "Vector::SetToLength(): length_after: " << Length() << std::endl;
      }
      else
      {
    	  std::cout << "==> nothing done in " << __FUNCTION__ << ": useless values: length: " << length << " value: " << MAX << std::endl;
#endif
      }
      return *this;
  }


  Vector &Vector::SetNonZeroElementsToMinusOrPlusOne()
  {
      for( std::vector< float>::iterator itr = this->begin(); itr != this->end(); ++itr)
      {
          if( *itr > 0)
          {
              *itr = 1;
          }
          else if( *itr < 0)
          {
              *itr = -1;
          }
      }
      return *this;
  }

  Vector &Vector::Merge( const Vector &V)
  {
      this->insert( this->end(), V.begin(), V.end());
      return *this;
  }

  bool Vector::AreAllElements( const float &VALUE) const
  {
      for( std::vector< float>::const_iterator itr = this->begin(); itr != this->end(); ++itr)
    	  if( *itr != VALUE)
		  {
    		  return false;
		  }
      return true;
  }

//  Vector &Vector::RemoveDuplicates()
//  {
//      std::vector< float> tmp;
//      std::vector< float>::iterator new_end = std::unique( this->begin(), this->end());
//      tmp.insert( tmp.end(), this->begin(), new_end);
//      *this = tmp;
//      return *this;
//  }

} // end namespace math

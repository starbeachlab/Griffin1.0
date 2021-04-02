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


#ifndef MATRIX_VECTOR_FUNCTIONS_H
#define MATRIX_VECTOR_FUNCTIONS_H

#include "matrix.h"

namespace math
{
  Matrix MatrixProduct( const std::vector< float> &V1, const std::vector< float> &V2);


  std::vector< float> Linearize( const Matrix &M);

  template< class ForwardIterator>
  void Linearize( const Matrix &M, ForwardIterator &ITR)
  {
	  for( size_t i = 0; i < M.GetNrRows(); ++i)
		  for( size_t j = 0; j < M.GetNrCols(); ++j)
		  {
			  *ITR++ = M( i, j);
		  }
  }



} // namespace math


#endif

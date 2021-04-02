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


#include "../../include/math/matrix_vector_functions.h"

namespace math
{
	Matrix MatrixProduct( const std::vector< float> &V1, const std::vector< float> &V2)
	{
		DebugWrite( __FUNCTION__);
		Matrix matrix( V1.size(), V2.size());
		for( size_t i( 0); i < V1.size(); ++i)
			for( size_t j( 0); j < V2.size(); ++j)
			{
				DebugWriteNoFlush( "(" << V1[i] << " * " << V2[j] << " = " << V1[i] * V2[j] << ") ");
				matrix( i, j) = V1[i] * V2[j];
			}
		DebugWrite( "");
		return matrix;
	}


	  std::vector< float> Linearize( const Matrix &M)
	  {
		  std::vector< float> v( M.GetNrCols() * M.GetNrRows());
		  std::vector< float>::iterator v_itr = v.begin();
		  for( size_t i = 0; i < M.GetNrRows(); ++i)
			  for( size_t j = 0; j < M.GetNrCols(); ++j)
			  {
				  *v_itr++ = M( i, j);
			  }
		  return v;
	  }



} // end namespace math

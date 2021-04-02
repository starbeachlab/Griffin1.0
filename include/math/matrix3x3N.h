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


#ifndef MATRIX3x3N_H
#define MATRIX3x3N_H

#include "matrix.h"

namespace math
{
	class Matrix3x3N
	: public Matrix
	{
	public:

		Matrix3x3N()
		: Matrix( 3, 3)
	    {}
      
		template< class ForwardIterator>
		Matrix3x3N( ForwardIterator &PTR)
		: Matrix( 3, 3)
		{
			Matrix::m_Data[0][0] = *PTR++;
			Matrix::m_Data[0][1] = *PTR++;
			Matrix::m_Data[0][2] = *PTR++;
			Matrix::m_Data[1][0] = *PTR++;
			Matrix::m_Data[1][1] = *PTR++;
			Matrix::m_Data[1][2] = *PTR++;
			Matrix::m_Data[2][0] = *PTR++;
			Matrix::m_Data[2][1] = *PTR++;
			Matrix::m_Data[2][2] = *PTR++;
		}

		virtual ~Matrix3x3N(){};

	}; // class Matrix3x3N

} // namespace math
#endif

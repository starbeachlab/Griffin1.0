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


#include "../../include/math/matrix.h"

namespace math
{

	float &Matrix::operator () ( const size_t &ROW, const size_t &COL)
	{
		assert( ROW < m_NrRows && COL < m_NrColumns);
		return m_Data[ ROW][ COL];
	}

	const float &Matrix::operator () ( const size_t &ROW, const size_t &COL) const
	{
		assert( ROW < m_NrRows && COL < m_NrColumns);
		return m_Data[ ROW][ COL];
	}

	Matrix &Matrix::operator += ( const Matrix &MATRIX)
	{
		assert( m_NrRows == MATRIX.m_NrRows && m_NrColumns == MATRIX.m_NrColumns);
		std::vector< std::vector< float> >::const_iterator m_itr = MATRIX.m_Data.begin();
		std::vector< float>::const_iterator m_jtr;
		for( std::vector< std::vector< float> >::iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr, ++m_itr)
		{
			m_jtr = m_itr->begin();
			for( std::vector< float>::iterator jtr = itr->begin(); jtr != itr->end(); ++jtr, ++m_jtr)
			{
				*jtr += *m_jtr;
			}
		}
		return *this;
	}


	Matrix &Matrix::operator -= ( const Matrix &MATRIX)
	{
		DebugWrite( __FUNCTION__);
		assert( m_NrRows == MATRIX.m_NrRows && m_NrColumns == MATRIX.m_NrColumns);
		std::vector< std::vector< float> >::const_iterator m_itr = MATRIX.m_Data.begin();
		std::vector< float>::const_iterator m_jtr;
		for( std::vector< std::vector< float> >::iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr, ++m_itr)
		{
			m_jtr = m_itr->begin();
			for( std::vector< float>::iterator jtr = itr->begin(); jtr != itr->end(); ++jtr, ++m_jtr)
			{
				DebugWriteNoFlush( "(" << *jtr << " - " << *m_jtr << " = ");
				*jtr -= *m_jtr;
				DebugWriteNoFlush( *jtr << ") ");
			}
		}
		DebugWrite("");
		return *this;
	}



} // end namespace math

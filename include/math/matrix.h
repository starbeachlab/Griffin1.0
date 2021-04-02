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


#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"

namespace math
{
  class Matrix
  {
  private:
    size_t  m_NrRows;
    size_t  m_NrColumns;
  protected:
    std::vector< std::vector< float> > m_Data;

  public:
    Matrix()
      : m_NrRows(0),
      m_NrColumns( 0),
      m_Data( 0, std::vector< float>( 0))
    {}

      Matrix( const size_t NR_ROWS, const size_t NR_COLS)
    : m_NrRows( NR_ROWS),
    m_NrColumns( NR_COLS),
	m_Data( NR_ROWS, std::vector< float>( NR_COLS))
      {}
    
    Matrix( const Matrix &MATRIX)
      : m_NrRows( MATRIX.m_NrRows),
      m_NrColumns( MATRIX.m_NrColumns),
      m_Data( MATRIX.m_Data)
        {}
      
      virtual ~Matrix(){}
      
      virtual float &operator () ( const size_t &ROW, const size_t &COL);
      virtual const float &operator () ( const size_t &ROW, const size_t &COL) const;
        
      virtual Matrix &operator += ( const Matrix &MATRIX);

      virtual Matrix &operator -= ( const Matrix &MATRIX);

      const size_t &GetNrCols() const{ return m_NrColumns;}

      const size_t &GetNrRows() const{ return m_NrRows;}


//
//      //! write matrix to ostream
//      virtual std::ostream& Write( std::ostream &STREAM) const;
//
//      //! read matrix from istream
//      virtual std::istream& Read( std::istream &STREAM);


  }; // class Matrix

//  inline std::ostream & operator << ( std::ostream &STREAM, const Matrix &M)
//  {
//      return M.Write( STREAM);
//  }
                      
} // namespace math
#endif

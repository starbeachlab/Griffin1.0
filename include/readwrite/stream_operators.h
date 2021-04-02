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


#ifndef STREAM_OPERATORS_H
#define STREAM_OPERATORS_H

#include <fstream>

class StreamOperators
{
 public:
  virtual ~StreamOperators(){}
  
  virtual std::ostream &Write( std::ostream &STREAM) const = 0;

  virtual std::istream &Read( std::istream &STREAM) = 0;
};

inline
 std::ostream& operator << ( std::ostream& STREAM, const StreamOperators &OPERATOR)
   { return OPERATOR.Write( STREAM);}
 
inline
 std::istream& operator >> ( std::istream& STREAM, StreamOperators &OPERATOR)
   { return OPERATOR.Read( STREAM);}

#endif
 

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


class SurfaceReader
{
public:
  

  SurfaceContainer
  ReadSurfaceFile( std::istream &STREAM)
  {
    std::string str;
    std::getline( STREAM, str);
    std::getline( STREAM, str);

    size_t nr_of_lines;
    
    STREAM >> nr_of_lines >> str >> str >> str;

    SurfaceContainer surface;

    for( size_t i = 0; i < nr_of_lines; ++i)
      {
    float nx, ny, nz, x, y, z;
    STREAM >> x >> y >> z >> nx >> ny >> nz >> str >> str >> str;

    surface.push_back( SurfacePtr( new SurfacePoint( x, y, z, nx, ny, nz)));
      }
    return surface;
  }
  

};  // end class LipidAtomReader

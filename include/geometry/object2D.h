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


#ifndef GEOM_OBJECT2D_H
#define GEOM_OBJECT2D_H

#include <vector>

#include "object.h"


namespace geom
{

  class Object2D
    : public Object
  {
  protected:
    math::Vector3N  m_XAxis;
    math::Vector3N  m_YAxis;
    
  public:
    
    Object2D()
      : Object(),
      m_XAxis(),
      m_YAxis()
      {}

    Object2D( const math::Vector3N &POS, const math::Vector3N &XAXIS, const math::Vector3N &YAXIS)
      : Object( POS),
      m_XAxis( XAXIS),
      m_YAxis( YAXIS)
      {}

    Object2D( const Object2D &OBJECT)
      : Object( OBJECT),
      m_XAxis( OBJECT.m_XAxis),
      m_YAxis( OBJECT.m_YAxis)
      {}


    virtual ~Object2D(){}
    virtual Object2D *Clone() const = 0;
    virtual math::Vector3N GetPosition() const = 0;
    virtual bool::IsPointWithinObject( const math::Vector3N) const;
    virtual bool::IsPointOnSurface( const math::Vector3N) const = 0;
    virtual std::vector< math::Vector3N> SurfacePointCloud( const float &RESOLUTION) const = 0;
    virtual std::vector< float> GetSizes() const;
  }; // end class Object2D

} // end namespace geom

#endif

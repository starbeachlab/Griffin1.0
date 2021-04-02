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


#ifndef GEOM_Object1D_H
#define GEOM_Object1D_H

#include <vector>

#include "object.h"


namespace geom
{

  class Object1D
    : public Object
  {
  protected:
    math::Vector3N  m_ZAxis;

  public:

    Object1D()
      : Object(),
      m_ZAxis()
      {}

    Object1D( const math::Vector3N &POS, const math::Vector3N &Z_AXIS = math::Vector3N( 0.0, 0.0, 1.0))
      : Object( POS),
      m_ZAxis( Z_AXIS.NormalizedCopy())
      {}

    Object1D( const Object1D &OBJECT)
      : Object( OBJECT),
      m_ZAxis( OBJECT.m_ZAxis)
      {}


    virtual ~Object1D(){}
    virtual Object1D *Clone() const = 0;
//    virtual math::Vector3N GetPosition() const = 0;

    virtual bool IsPointWithin( const math::Vector3N) const
    {
        return false;
    }

    virtual std::vector< math::Vector3N> GetAxes() const
    {
        return std::vector< math::Vector3N> ( 1, m_ZAxis);
    }


//    virtual bool IsPointOnSurface( const math::Vector3N) const = 0;
//    virtual std::vector< math::Vector3N> SurfacePointCloud( const float &RESOLUTION) const = 0;
//    virtual std::vector< float> GetSizes() const;
  }; // end class Object1D

} // end namespace geom

#endif

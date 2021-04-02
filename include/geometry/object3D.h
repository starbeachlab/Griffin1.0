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


#ifndef GEOM_OBJECT3D_H
#define GEOM_OBJECT3D_H

#include <vector>

#include "object.h"


namespace geom
{

  class Object3D
    : public Object
  {
  protected:
    math::Vector3N  m_XAxis;
    math::Vector3N  m_YAxis;
    math::Vector3N  m_ZAxis;

  public:

    Object3D()
      : Object( math::Vector3N( 0.0, 0.0, 0.0)),
      m_XAxis( 1.0, 0.0, 0.0),
      m_YAxis( 0.0, 1.0, 0.0),
      m_ZAxis( 0.0, 0.0, 1.0)
      {}

    Object3D( const math::Vector3N &POS)
      : Object( POS),
      m_XAxis( 1.0, 0.0, 0.0),
      m_YAxis( 0.0, 1.0, 0.0),
      m_ZAxis( 0.0, 0.0, 1.0)
      {}

    Object3D( const math::Vector3N &POS, const math::Vector3N &XAXIS, const math::Vector3N &YAXIS, const math::Vector3N &ZAXIS)
      : Object( POS),
      m_XAxis( XAXIS),
      m_YAxis( YAXIS),
      m_ZAxis( ZAXIS)
      {}

    Object3D( const Object3D &OBJECT)
      : Object( OBJECT),
      m_XAxis( OBJECT.m_XAxis),
      m_YAxis( OBJECT.m_YAxis),
      m_ZAxis( OBJECT.m_ZAxis)
      {}




    virtual ~Object3D(){}
    virtual Object3D *Clone() const = 0;


    virtual std::vector< math::Vector3N> GetAxes() const
    {
        std::vector< math::Vector3N> v( 3);
        v[0] = m_XAxis;
        v[1] = m_YAxis;
        v[2] = m_ZAxis;
        return v;
    }


    virtual math::Vector3N GetAxes( const size_t &ID) const
    {
        switch( ID)
        {
        case 0: return m_XAxis;
        case 1: return m_YAxis;
        case 2: return m_ZAxis;
        default: return math::Vector3N();
        }
    }


    virtual bool IsPointWithin( const math::Vector3N &POS) const{ return false;}

  }; // end class Object3D

} // end namespace geom

#endif

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


#ifndef PLANE_FUNCTIONS_H_
#define PLANE_FUNCTIONS_H_

#include "plane.h"

namespace geom
{

    inline float DistancePointToPlane( const math::Vector3N &POSITION, const Plane &PLANE)
    {
        math::Vector3N
            connection_point_origin( POSITION - PLANE.GetPosition());
        return std::abs( math::ScalarProduct( connection_point_origin, PLANE.GetNormalVector()) / PLANE.GetNormalVector().Length());
    }

    inline std::pair< math::Vector3N, float> FootpointAndDistanceOfPointToPlane( const math::Vector3N &POSITION, const Plane &PLANE)
    {
        math::Vector3N
            norm( PLANE.GetNormalVector()),
            connection_point_origin( POSITION - PLANE.GetPosition());
        float
            angle( math::Angle( norm, connection_point_origin));

        if( angle > math::HalfPi)
        {
            norm *= -1.0;
            angle = math::Pi - angle;
        }

        float
            dist( connection_point_origin.Length() * cos( angle));
        math::Vector3N
            footpoint( POSITION - dist * norm);
        return std::make_pair( footpoint, dist);
    }

} // end namespace geom




#endif /* PLANE_FUNCTIONS_H_ */

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


#ifndef TORUS_FUNCTIONS_H_
#define TORUS_FUNCTIONS_H_

#include "torus_slice.h"
#include "circle_functions.h"

namespace geom
{
    inline
    bool
    IsPointSurroundedByTorusSlice( const math::Vector3N &POSITION, const TorusSlice &TORUS, const float &THRESHOLD = 0.0)
    {
    //        std::cout << __PRETTY_FUNCTION__ << std::endl;
        math::Vector3N
            connection_point_origin( POSITION - TORUS.GetPosition());
        float
            distance_point_origin( connection_point_origin.Length());

        if( distance_point_origin > TORUS.GetRadius())
        {
            return false;
        }

        float
            pos_or_neg_distance_point_to_plane( math::ScalarProduct( connection_point_origin, TORUS.GetNormalVector())),
            distance_minimum_angle_to_plane( TORUS.GetTorusRadius() * sin( TORUS.GetMinimumAngle())),
            distance_maximum_angle_to_plane( TORUS.GetTorusRadius() * sin( TORUS.GetMaximumAngle()));

        bool
            is_within_slice
            (
                    ( pos_or_neg_distance_point_to_plane > distance_minimum_angle_to_plane && pos_or_neg_distance_point_to_plane < distance_maximum_angle_to_plane)
                    || ( pos_or_neg_distance_point_to_plane < distance_minimum_angle_to_plane && pos_or_neg_distance_point_to_plane > distance_maximum_angle_to_plane)
            );

        if( !is_within_slice)
        {
            return false;
        }

        float
            distance_point_plane( DistancePointToPlane( POSITION, TORUS)),
            distance_projection_origin( std::sqrt( math::Square( distance_point_origin) - math::Square( distance_point_plane))),
            difference_distance_radius( distance_projection_origin - TORUS.GetRadius()),
            distance_point_circle( std::sqrt( math::Square( distance_point_plane) + math::Square( difference_distance_radius)));

        //        return ( distance_point_circle - TORUS.GetTorusRadius()r < 0 && distance_point_circle - TORUS.GetTorusRadius() > -0.2); // surface points
            return ( distance_point_circle  > TORUS.GetTorusRadius() + THRESHOLD);
    }

    inline
    bool
    IsPointSurroundedByTorus( const math::Vector3N &POSITION, const Torus &TORUS, const float &THRESHOLD = 0.0)
    {
    //        std::cout << __PRETTY_FUNCTION__ << std::endl;
        float
            distance_point_origin( ( POSITION - TORUS.GetPosition()).Length()),
            distance_point_plane( DistancePointToPlane( POSITION, TORUS)),
            distance_projection_origin( std::sqrt( math::Square( distance_point_origin) - math::Square( distance_point_plane))),
            difference_distance_radius( distance_projection_origin - TORUS.GetRadius()),
            distance_point_circle( std::sqrt( math::Square( distance_point_plane) + math::Square( difference_distance_radius)));

    //        return ( distance_point_circle - TORUS.GetTorusRadius()r < 0 && distance_point_circle - TORUS.GetTorusRadius() > -0.2); // surface points
        return ( distance_point_circle  > TORUS.GetTorusRadius() + THRESHOLD && distance_point_plane <= TORUS.GetTorusRadius() && distance_point_origin < TORUS.GetRadius());
    }

    inline
    bool
    IsPointOnTorusSurfaceWithinThreshold(  const math::Vector3N &POSITION, const Torus &TORUS, const float &PRECISION)
    {
        float
            distance_point_circle( DistancePointToCircle( POSITION, TORUS));

        return ( distance_point_circle - TORUS.GetTorusRadius() < 1e-5 && distance_point_circle - TORUS.GetTorusRadius() > -PRECISION); // surface points
    }

    inline
    bool
    IsPointWithinTorus( const math::Vector3N &POSITION, const Torus &TORUS)
    {
        return ( DistancePointToCircle( POSITION, TORUS) < TORUS.GetTorusRadius() + 1e-5);
    }

    inline
    float
    DistanceBetweenTwoTori( const Torus &A, const Torus &B)
    {
        return DistanceBetweenTwoCircles( A, B) - A.GetTorusRadius() - B.GetTorusRadius();
    }

    inline
    bool
    DoToriOverlap( const Torus &A, const Torus &B)
    {
        return DistanceBetweenTwoTori( A, B) == 0;
    }

} // end namespace geom




#endif /* TORUS_FUNCTIONS_H_ */

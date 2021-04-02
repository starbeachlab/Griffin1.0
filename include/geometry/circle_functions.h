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


#ifndef CIRCLE_FUNCTIONS_H_
#define CIRCLE_FUNCTIONS_H_

#include "circle.h"
#include "plane_functions.h"

namespace geom
{

    inline
    float
    DistancePointToCircle( const math::Vector3N &POSITION, const Circle &CIRCLE)
    {
//        std::cout << __PRETTY_FUNCTION__ << std::endl;
        float
            distance_point_origin( ( POSITION - CIRCLE.GetPosition()).Length()),
            distance_point_plane( DistancePointToPlane( POSITION, CIRCLE)),
            distance_projection_origin( std::sqrt( math::Square( distance_point_origin) - math::Square( distance_point_plane))),
            difference_distance_radius( distance_projection_origin - CIRCLE.GetRadius());

        return std::sqrt( math::Square( distance_point_plane) + math::Square( difference_distance_radius));
    }

    inline
    float
    DistanceBetweenTwoCircles( const Circle &A, const Circle &B)
    {
        // either go for old line segment implementation
        math::Vector3N
            connection( B.GetPosition() - A.GetPosition()),
            a_y( math::CrossProduct( A.GetNormalVector(), connection).Normalize()),
            a_x( math::CrossProduct( a_y, A.GetNormalVector()).Normalize()),
            b_y( math::CrossProduct( B.GetNormalVector(), -1.0 * connection).Normalize()),
            b_x( math::CrossProduct( b_y, B.GetNormalVector()).Normalize());

        float
            volume( math::TripleProduct( connection, A.GetNormalVector(), a_x));
        if( volume > 1e-5)
        {
            std::cout << "===> too much volume for first: " << volume << std::endl;
        }
        volume = math::TripleProduct( connection, B.GetNormalVector(), b_x);
        if( volume > 1e-5)
        {
            std::cout << "===> too much volume for second: " << volume << std::endl;
        }

        float
            angle( math::RadiansToDegrees( math::Angle( a_x, connection)));
        if( angle > 90)
        {
            std::cout << "===> angle first: " << angle << std::endl;
        }
        angle = math::RadiansToDegrees( math::Angle( b_x, connection));
        if( angle < 90)
        {
            std::cout << "===> angle second: " << angle << std::endl;
        }

//        std::cout << "size first: " << a_x.Length() << std::endl;
//        std::cout << "size second: " << b_x.Length() << std::endl;

        // or go for the laws of sinus and cosinus (might be faster)
        float
            alpha_test( math::Angle( a_x, b_x));

        float
            distance( connection.Length()),
            a_angle( math::Angle( connection, a_x)), // angle opposite of triangle side belonging to torus A
            b_angle( math::Angle( -1.0 * connection, b_x)), // angle opposite of triangle side belonging to torus B
            distance_triangular;

        if( a_angle > math::HalfPi)
        {
            a_angle = math::Pi - a_angle;
        }
        if( b_angle > math::HalfPi)
        {
            b_angle = math::Pi - b_angle;
        }
        if( fabs( math::Pi - a_angle - b_angle - alpha_test) > 1e-2)
        {
            std::cout << "===> angle mismatch: pi - " << math::RadiansToDegrees( a_angle) << " - " << math::RadiansToDegrees( b_angle) << " = " << math::RadiansToDegrees( math::Pi - a_angle - b_angle) << " != " << math::RadiansToDegrees( alpha_test) << std::endl;
        }

        float
            a_side_length( math::LengthFromLawOfSinus( a_angle, distance, alpha_test)), // triangle side belonging to torus A
            b_side_length( math::LengthFromLawOfSinus( b_angle, distance, alpha_test)); // triangle side belonging to torus B

        static size_t ccc( 0);
        std::cout << "mol new" << std::endl;
        std::cout << "graphics " << ccc << " color white" << std::endl;
        std::cout << "graphics " << ccc << " cone {" <<  A.GetPosition()( 0)<< " " << A.GetPosition()( 1) << " " << A.GetPosition()( 2) << "} {" << A.GetPosition()( 0) + A.GetNormalVector()( 0) << "  " << A.GetPosition()( 1) + A.GetNormalVector()( 1) << "  " << A.GetPosition()( 2) + A.GetNormalVector()( 2) << "} radius " << A.GetRadius() << " resolution 30 " << std::endl;
        std::cout << "graphics " << ccc << " cone {" <<  B.GetPosition()( 0)<< " " << B.GetPosition()( 1) << " " << B.GetPosition()( 2) << "} {" << B.GetPosition()( 0) + B.GetNormalVector()( 0) << "  " << B.GetPosition()( 1) + B.GetNormalVector()( 1) << "  " << B.GetPosition()( 2) + B.GetNormalVector()( 2) << "} radius " << B.GetRadius()  << " resolution 30 " << std::endl;
        std::cout << "graphics " << ccc << " color blue" << std::endl;
        std::cout << "graphics " << ccc << " cylinder {" << A.GetPosition()( 0)    << " " << A.GetPosition()( 1) << " " << A.GetPosition()( 2) << "} {" << B.GetPosition()( 0) <<  "  " << B.GetPosition()( 1) << "  " << B.GetPosition()( 2) << "} radius 0.05" << std::endl;

        std::pair< math::Vector3N, float> pair;

        if( a_side_length < A.GetRadius() && b_side_length < B.GetRadius())
        {
            distance_triangular = 0.0;
        }
        else if( a_side_length < A.GetRadius())
        {
            pair = FootpointAndDistanceOfPointToPlane( B.GetPosition() + B.GetRadius() * b_x, A);
            distance_triangular = pair.second;

//            math::Vector3N pos( B.GetPosition() + B.GetRadius() * b_x), end( pair.first);
            math::Vector3N pos( A.GetPosition() + A.GetRadius() * a_x), end( pair.first);
            std::cout << "graphics " << ccc << " color orange" << std::endl;
            std::cout << "graphics " << ccc << " cone {" << pos( 0)    << " " << pos( 1) << " " << pos( 2) << "} {" << end( 0) <<  "  " << end( 1) << "  " << end( 2) << "} radius 0.05" << std::endl;
        }
        else if( b_side_length < B.GetRadius())
        {
            pair = FootpointAndDistanceOfPointToPlane( A.GetPosition() + A.GetRadius() * a_x, B);
            distance_triangular = pair.second;

//            math::Vector3N pos( A.GetPosition() + A.GetRadius() * a_x), end( pair.first);
            math::Vector3N pos( B.GetPosition() + B.GetRadius() * b_x), end( pair.first);
            std::cout << "graphics " << ccc << " color green" << std::endl;
            std::cout << "graphics " << ccc << " cone {" << pos( 0)    << " " << pos( 1) << " " << pos( 2) << "} {" << end( 0) <<  "  " << end( 1) << "  " << end( 2) << "} radius 0.05" << std::endl;

        }
        else
        {
            distance_triangular = math::Distance( A.GetPosition() + A.GetRadius() * a_x, B.GetPosition() + B.GetRadius() * b_x);

            math::Vector3N pos( A.GetPosition() + A.GetRadius() * a_x), end( B.GetPosition() + B.GetRadius() * b_x);
            if( fabs( distance_triangular - ( pos-end).Length()) > 1e-5)
            {
                std::cout << "===> mismatch distances: " << distance_triangular << " vs. " << (pos-end).Length() << std::endl;
            }
            std::cout << "graphics " << ccc << " color red" << std::endl;
            std::cout << "graphics " << ccc << " cylinder {" << pos( 0)    << " " << pos( 1) << " " << pos( 2) << "} {" << end( 0) <<  "  " << end( 1) << "  " << end( 2) << "} radius 0.05" << std::endl;
        }
        angle = math::RadiansToDegrees( math::Angle( a_x, connection));
        if( angle > 90)
        {
            std::cout << "===> angle first: " << angle << std::endl;
        }
        angle = math::RadiansToDegrees( math::Angle( b_x, connection));
        if( angle > 90)
        {
            std::cout << "===> angle second: " << angle << std::endl;
        }
        float min( distance_triangular);
        math::Vector3N min_a, min_b;
        std::cout << "graphics " << ccc << " color yellow" << std::endl;
        for( float x = 0; x < math::TwoPi; x += math::DegreesToRadians( 5.0))
            for( float y = 0; y < math::TwoPi; y+= math::DegreesToRadians( 5.0))
            {
                math::Vector3N pos_a( A.GetPosition() + A.GetRadius() * ( cos( x) * a_x + sin( x) * a_y));
                math::Vector3N pos_b( B.GetPosition() + B.GetRadius() * ( cos( y) * b_x + sin( y) * b_y));
                float disty( math::Distance( pos_a, pos_b));
                if( disty < min)
                {
                    min = disty;
                    min_a = pos_a;
                    min_b = pos_b;
                }
            }
        std::cout << "graphics " << ccc << " cylinder {" << min_a( 0)    << " " << min_a( 1) << " " << min_a( 2) << "} {" << min_b( 0) <<  "  " << min_b( 1) << "  " << min_b( 2) << "} radius 0.02" << std::endl;

//        if( fabs( distance_line_segment - distance_triangular) > 1e-2)
//        {
//            std::cout << "===> distances don't match: " << distance_line_segment << " " << distance_triangular << std::endl;
//        }


        ++ccc;

        return distance_triangular;
    }

//    inline
//    float
//    DistanceBetweenTwoCircles( const Circle &A, const Circle &B)
//    {
//        // initial guess
//        // if b is on top of a - solved
//        // else: force for direction
//        // minimize
//    }
} // end namespace geom




#endif /* CIRCLE_FUNCTIONS_H_ */

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


#ifndef POINT_SURFACE_CYLINDER_H_
#define POINT_SURFACE_CYLINDER_H_


#include "cylinder.h"
#include "directed_point_surface.h"

namespace geom
{

    class PointSurfaceCylinder
    : public DirectedPointSurfaceObject
    {

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        PointSurfaceCylinder()
        : DirectedPointSurfaceObject()
        {}

        //! construct from data
        PointSurfaceCylinder( const math::Vector3N &POS, const math::Vector3N &ZAXIS, const float &RADIUS, const float &HEIGHT, const float &RESOLUTION)
        : DirectedPointSurfaceObject( boost::shared_ptr< Object>( new Cylinder( POS, ZAXIS, RADIUS, HEIGHT)))
        { BuildSurface( RESOLUTION);}

        //! construct from data
        PointSurfaceCylinder( const math::Vector3N &POS, const float &RADIUS, const float &HEIGHT, const float &RESOLUTION)
        : DirectedPointSurfaceObject( boost::shared_ptr< Object>( new Cylinder( POS, RADIUS, HEIGHT)))
        { BuildSurface( RESOLUTION);}

        //! copy constructor
        PointSurfaceCylinder( const PointSurfaceCylinder &ORIGINAL)
        : DirectedPointSurfaceObject( ORIGINAL)
        {}

        //! virtual destructor
        virtual ~PointSurfaceCylinder(){}

        //! virtual copy constructor
        virtual PointSurfaceCylinder *Clone() const{ return new PointSurfaceCylinder( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        float RingSectionSurface( const float &INNER_RADIUS, const float &OUTER_RADIUS, const float &ANGLE)
        {
            return ANGLE * 0.5 * ( INNER_RADIUS + OUTER_RADIUS) * ( OUTER_RADIUS - INNER_RADIUS);
        }


        float AngleFromRingSectionSurface( const float &INNER_RADIUS, const float &OUTER_RADIUS, const float &SURFACE)
        {
            return 2.0 * SURFACE / ( ( INNER_RADIUS + OUTER_RADIUS) * ( OUTER_RADIUS - INNER_RADIUS));
        }

        float DetermineDeltaToBeCloseToGuess( const float &INTERVAL, const float &INITIAL_DELTA)
        {
            return INTERVAL / ( float( 1 + size_t( INTERVAL / INITIAL_DELTA)));
        }

    protected:
        void BuildSurface( const float &RESOLUTION)
        {
            if( RESOLUTION > math::TwoPi)
            {
                std::cout << "===> check resolution in " << GetClassName() << " " << __FUNCTION__ << ": " << RESOLUTION << " (max: " << math::TwoPi << ")" << std::endl;
            }
            std::vector< float> sizes( GetSizes());

            float
                radius( sizes[0]),
                height( sizes[1]),
                dist( RESOLUTION * radius),
                delta_z( DetermineDeltaToBeCloseToGuess( height, dist)),
                delta_r( DetermineDeltaToBeCloseToGuess( radius, dist)),
                standard_surf( RingSectionSurface( radius - delta_r, radius, RESOLUTION)),
                r,
                z,
                angle,
                delta_angle = -8888.88;

            math::Vector3N
                z_axis( GetAxes()[ 0].Normalize()),
                y_axis( math::CrossProduct( z_axis, math::Vector3N( 1.0, 0.0, 0.0))),
                x_axis( math::CrossProduct( y_axis, z_axis).Normalize()),
                position,
                norm;

            y_axis = math::CrossProduct( z_axis, x_axis).Normalize();


            // upper & lower circle
            for( int i = -1; i <= 1; i += 2)
            {
                z = float( i) * 0.5 * height;
                // top and bottom
                for( r = 0.0; r < radius; r += delta_r)
                {
                    dist = r + 0.5 * delta_r;
                    delta_angle = AngleFromRingSectionSurface( r, r + delta_r, standard_surf);
                    delta_angle =  DetermineDeltaToBeCloseToGuess( math::TwoPi, delta_angle);
                    for( angle = 0.0; angle <= math::TwoPi; angle += delta_angle)
                    {
                        position = GetPosition() + z * z_axis + dist * ( cos( angle) * x_axis + sin( angle) * y_axis);
                        PushBackDirectedSurfacePoint( boost::shared_ptr< DirectedSurfacePoint>( new DirectedSurfacePoint( position, float( i) * z_axis)));
                    }
                }
                // edges
                for( angle = 0.0; angle <= math::TwoPi; angle += delta_angle)
                {
                    norm = cos( angle) * x_axis + sin( angle) * y_axis;
                    position = GetPosition() + z * z_axis + radius * norm;
                    PushBackDirectedSurfacePoint( boost::shared_ptr< DirectedSurfacePoint>( new DirectedSurfacePoint( position, float( i) * z_axis + norm)));
                }
            }

            // mantel
            for( z = -0.5 * height + delta_z; z < 0.5 * height; z += delta_z)
            {
                for( angle = 0.0; angle <= math::TwoPi; angle += delta_angle)
                {
                    norm = cos( angle) * x_axis + sin( angle) * y_axis;
                    position = GetPosition() + z * z_axis + radius * norm;
                    PushBackDirectedSurfacePoint( boost::shared_ptr< DirectedSurfacePoint>( new DirectedSurfacePoint( position, norm)));
                }
            }
        }

        /////////////////////////
        //      Read/Write     //
        /////////////////////////
    public:
//        virtual std::istream &Read( std::istream &STREAM)
//        {
//            return PointSurfaceObject::Read( STREAM);
//        }
//
//        virtual std::ostream &Write( std::ostream &STREAM) const
//        {
//            return PointSurfaceObject::Write( STREAM);
//        }
//        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM) const
//        { return PointSurfaceObject::WriteVmdCommands( STREAM);}

//    //        STREAM << "mol new" << std::endl;
//            STREAM << "graphics 0 color red" <<  std::endl;
//            STREAM << "graphics 0 sphere {" << GetPosition()(0) << " " << GetPosition()(1) << " " << GetPosition()(2) << "} radius 0.1 resolution 21" << std::endl;
//            return STREAM;
//        }

        virtual std::string GetClassName() const
        { return mystr::GetClassName( __PRETTY_FUNCTION__);}

    }; // end class PointSurfaceCylinder
} // end namespace geom




#endif /* POINT_SURFACE_CYLINDER_H_ */

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


#ifndef POINT_SURFACE_CUBE_H_
#define POINT_SURFACE_CUBE_H_


#include "cube.h"
#include "point_surface.h"
//#include ""

namespace geom
{

    class PointSurfaceCube
    : public DirectedPointSurfaceObject
    {

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        PointSurfaceCube()
        : DirectedPointSurfaceObject()
        {}

        //! construct from data
        PointSurfaceCube( const math::Vector3N &POS, const float &XSIZE, const float &YSIZE, const float &ZSIZE, const float &RESOLUTION)
        : DirectedPointSurfaceObject( boost::shared_ptr< Object>( new Cube( POS, XSIZE, YSIZE, ZSIZE)))
        {
            DebugWrite( __FUNCTION__);
            BuildCubeSurface( RESOLUTION);
        }

        //! copy constructor
        PointSurfaceCube( const PointSurfaceCube &ORIGINAL)
        : DirectedPointSurfaceObject( ORIGINAL)
        {}

        //! virtual destructor
        virtual ~PointSurfaceCube(){}

        //! virtual copy constructor
        virtual PointSurfaceCube *Clone() const{ return new PointSurfaceCube( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

    protected:
        void BuildCubeSurface( const float &RESOLUTION)
        {
            DebugWrite( __FUNCTION__);
            math::Vector3N
                position,
                zero( 0.0, 0.0, 0.0),
                direction;

            for( size_t i = 0; i < 3; ++i)
            {
                math::Vector3N normal( zero);
                normal( i) = 1.0;
                size_t a( ( i + 1) % 3);
                size_t b( ( i + 2) % 3);
                std::vector< float> sizes( GetSizes());
                sizes *= 0.5;
                DebugWrite( "sizes size: " << sizes.size());
                for( float x = -1.0 * sizes[ a]; x < sizes[ a]; x += RESOLUTION)
                    for( float y = -1.0 * sizes[ b]; y < sizes[ b]; y += RESOLUTION)
                    {
                        position = DirectedPointSurfaceObject::GetPosition();
                        position( i) += sizes[ i];
                        position( a) += x;
                        position( b) += y;
                        direction = normal;
//                        DebugWrite( "direction: " << direction);
                        if( math::IsEqualWithinThreshold( x, -1.0 * sizes[ a], 0.01))
                        {
                            direction( a) -= 1.0;
//                            DebugWrite( "direction: " << direction);
                        }
                        if( math::IsEqualWithinThreshold( y, -1.0 * sizes[ b], 0.01))
                        {
                            direction( b) -= 1.0;
//                            DebugWrite( "direction: " << direction);
                        }
                        direction.Normalize();
                        float res( 0.01);  // NO PROPER SURFACE CALCULATION YET
                        PushBackDirectedSurfacePoint( boost::shared_ptr< DirectedSurfacePoint>( new DirectedSurfacePoint( position, direction, res)));
//                        position( i) = DirectedPointSurfaceObject::GetPosition()( i) - sizes[ i];
//                        direction( i) *= -1.0;
                        position( i) += -2.0 * sizes[ i];
                        position( a) += -2.0 * x;
                        position( b) += -2.0 * y;
                        direction *= -1.0;
                        PushBackDirectedSurfacePoint( boost::shared_ptr< DirectedSurfacePoint>( new DirectedSurfacePoint( position, direction, res)));
                    }
            }
        }

        //                for( size_t i = 0; i < 3; ++i)
//                {
//                    math::Vector3N normal( zero);
//                    normal( i) = 1.0;
//                    size_t a( ( i + 1) % 3);
//                    size_t b( ( i + 2) % 3);
//                    std::vector< float> sizes( GetSizes());
//                    sizes *= 0.5;
//                    DebugWrite( "sizes size: " << sizes.size());
//                for( float x = -1.0 * sizes[ a] + 0.5 * RESOLUTION; x < sizes[ a]; x += RESOLUTION)
//                    for( float y = -1.0 * sizes[ b] + 0.5 * RESOLUTION; y < sizes[ b]; y += RESOLUTION)
//                    {
//                        position = DirectedPointSurfaceObject::GetPosition();
//                        position( i) += sizes[ i];
//                        position( a) += x;
//                        position( b) += y;
//                        float res( 0.01);
//                        PushBackDirectedSurfacePoint( boost::shared_ptr< DirectedSurfacePoint>( new DirectedSurfacePoint( position, normal, res)));
//                        position( i) = DirectedPointSurfaceObject::GetPosition()( i) - sizes[ i];
//                        PushBackDirectedSurfacePoint( boost::shared_ptr< DirectedSurfacePoint>( new DirectedSurfacePoint( position, -1.0 * normal, 0.01)));  // NO PROPER SURFACE CALCULATION YET
//                    }

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

    }; // end class PointSurfaceCube
} // end namespace geom




#endif /* POINT_SURFACE_CUBE_H_ */

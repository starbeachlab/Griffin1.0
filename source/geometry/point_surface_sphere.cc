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


#include "../../include/geometry/point_surface_sphere.h"

namespace geom
{
    boost::scoped_ptr< DirectedPointSurface> PointSurfaceSphere::s_DefaultSurf( new DirectedPointSurface());

    //! default constructor
    PointSurfaceSphere::PointSurfaceSphere()
    : DirectedPointSurfaceObject()
    {
//        std::cout << mystr::GetClassName( __PRETTY_FUNCTION__) << "default constructor" << std::endl;
    }

    //! construct from position and radius
    PointSurfaceSphere::PointSurfaceSphere( const math::Vector3N &POS, const float &RADIUS, const float &RESOLUTION)
    : DirectedPointSurfaceObject( boost::shared_ptr< Object>( new Sphere( POS, RADIUS)), DirectedPointSurface( ( *s_DefaultSurf * RADIUS) + POS))
    {
        DebugWrite( __FUNCTION__ << " radius: " << RADIUS << " resolution: " << RESOLUTION);
//        std::cout << mystr::GetClassName( __PRETTY_FUNCTION__) << "construct from data" << std::endl;
        if( s_DefaultSurf->GetNrPoints() == 0)
        {
            DebugWrite( "initialize default sphere");
            InitializeDefaultSphere( RESOLUTION);
            DebugWrite( "initialize done");
            SetSurface( ( *s_DefaultSurf * GetRadius()) + POS);
            DebugWrite( "default surf set");

//            for( std::vector< boost::shared_ptr< SurfacePoint> >::const_iterator itr( s_DefaultSurf->GetData().begin()); itr != s_DefaultSurf->GetData().end(); ++itr)
//            {
////                std::cout << "check: ";
//                ( *itr)->Clone();
//            }
//            std::cout << "nr of points: " << tmp.GetNrPoints() << std::endl;
//            m_Surf = boost::shared_ptr< PointSurface>( new PointSurface( tmp));
        }
    }

    //! copy constructor
    PointSurfaceSphere::PointSurfaceSphere( const PointSurfaceSphere &ORIGINAL)
    : DirectedPointSurfaceObject( ORIGINAL)
    {
//        std::cout << mystr::GetClassName( __PRETTY_FUNCTION__) << " copy constructor" << std::endl;
    }

    //! virtual destructor
    PointSurfaceSphere::~PointSurfaceSphere(){}

    //! virtual copy constructor
    //!(needed for copying pointers to derived classes)
    PointSurfaceSphere *PointSurfaceSphere::Clone() const{ return new PointSurfaceSphere( *this);}


    const boost::scoped_ptr< DirectedPointSurface> &PointSurfaceSphere::InitializeDefaultSphere( const float &DELTA)
    {
        DebugWrite( __FUNCTION__);
//        std::cout << "initialize default sphere" << std::endl;
        int
            theta_max = 180,
            rings( theta_max / int( DELTA)),  // Adjust this in case Delta shall be < 1.
            nr_intervalls;
//            count;
        float
            theta_1( 0),
            theta_2( 0),
//            s_0( 0),
            delta_phi( 0),
            phi( 0),
            surf( 0),
            delta_radians( DELTA * math::Pi / float( 180.0));

          math::Vector3N position;

          // the initial surface area
          float surf_0( delta_radians * sin( delta_radians));// * r^2, but r = 1

          for( int i = 0 ; i < rings ; ++i )
          {
              // calculating delta_theta and delta_phi such that the surface area is as close to surf_0 as possible.
              theta_1 =  i * delta_radians;
              theta_2 = (i+1) * delta_radians;
              delta_phi = surf_0 / ( cos( theta_1 ) - cos( theta_2 )); // /r^2, but r = 1
              nr_intervalls = int( 2.0 * math::Pi / delta_phi) + 2;
              delta_phi = 2.0 * math::Pi / float(nr_intervalls);

              // the actual surface area, assigned to a regarded integration point
              surf =  delta_phi * ( cos( theta_1) - cos( theta_2) );  // * r^2, but r = 1

              // testing if a certain point on the surface of a sphere is overlapped by other spheres.
              for( int j = 0 ; j < nr_intervalls ; j ++ )
              {
                  phi = j * delta_phi;
                  position(0) = sin( ( i+0.5) * delta_radians ) * cos( phi); // * r, but r = 1
                  position(1) = sin( ( i+0.5) * delta_radians ) * sin( phi); // * r, but r = 1
                  position(2) = cos( ( i+0.5) * delta_radians );             // * r, but r = 1
                  s_DefaultSurf->PushBack( boost::shared_ptr< DirectedSurfacePoint>( new DirectedSurfacePoint( position, position.NormalizedCopy(), surf)));
              }
          }
          return s_DefaultSurf;
    }


    std::vector< math::Vector3N> PointSurfaceSphere::GetSurfacePointPositions() const  // or by
    {
        std::vector< math::Vector3N> pos( GetDirectedSurfacePoints().size());
        std::vector< math::Vector3N>::iterator pos_itr( pos.begin());
        std::vector< boost::shared_ptr< DirectedSurfacePoint> >::const_iterator surf_itr( GetDirectedSurfacePoints().begin());
        for( ; surf_itr != GetDirectedSurfacePoints().end(); ++surf_itr, ++pos_itr)
        {
            *pos_itr = ( *surf_itr)->GetPosition();
            //            *pos_itr = Sphere::GetPosition() + ( *surf_itr)->GetPosition();
        }
        return pos;
    }



//    std::ostream& PointSurfaceSphere::Write( std::ostream &STREAM) const
//    {
//        STREAM << mystr::GetClassName( __PRETTY_FUNCTION__) << std::endl;
//        Sphere::Write( STREAM);
//        STREAM << m_Surf << std::endl;
//        return STREAM;
//    }
//
//    std::ostream &PointSurfaceSphere::WriteVmdCommands( std::ostream &STREAM) const
//    {
////        STREAM << "mol new" << std::endl;
//        STREAM << "graphics 0 color red" <<  std::endl;
//        STREAM << "graphics 0 sphere {" << Sphere::GetPosition()(0) << " " << Sphere::GetPosition()(1) << " " << Sphere::GetPosition()(2) << "} radius 0.1 resolution 21" << std::endl;
//        m_Surf->WriteVmdCommands( STREAM);
//        return STREAM;
//    }
//
//
//    std::istream& PointSurfaceSphere::Read( std::istream &STREAM)
//    {
//        std::string str;
//        STREAM >> str;
//        assert( str == mystr::GetClassName( __PRETTY_FUNCTION__));
//        Sphere::Read( STREAM);
//        STREAM >> m_Surf;
//        return STREAM;
//    }

    std::string PointSurfaceSphere::GetClassName() const
    { return mystr::GetClassName( __PRETTY_FUNCTION__);}


    std::ostream &operator << ( std::ostream &STREAM, const PointSurfaceSphere &SPHERE)
    { return SPHERE.Write( STREAM);}

    std::istream &operator >> ( std::istream &STREAM, PointSurfaceSphere &SPHERE)
    { return SPHERE.Read( STREAM);}


} // end namespace geom




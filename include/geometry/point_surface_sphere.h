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


#ifndef POINT_SURFACE_SPHERE_H_
#define POINT_SURFACE_SPHERE_H_


#include "sphere.h"
#include "point_surface.h"
#include "../math/constants.h"
#include "directed_surface_point.h"
#include "../string/io_string_functions.h"
#include "../external/boost_functions.h"
#include "directed_point_surface_object.h"

namespace geom
{

    class PointSurfaceSphere
    : public DirectedPointSurfaceObject
    {
    protected:
    ///////////////
    //    DATA     //
    ///////////////
        static boost::scoped_ptr< DirectedPointSurface>  s_DefaultSurf;   //!< common surface of given resolution for radius 1
//        boost::shared_ptr< PointSurface>          m_Surf;


    public:
    //////////////////////////////////
    //  CONSTRUCTION & DESTRUCTION  //
    //////////////////////////////////

        //! default constructor
        PointSurfaceSphere();

        //! construct from position and radius
        PointSurfaceSphere( const math::Vector3N &POS, const float &RADIUS, const float &RESOLUTION);

        //! copy constructor
        PointSurfaceSphere( const PointSurfaceSphere &ORIGINAL);

        //! virtual destructor
        virtual ~PointSurfaceSphere();

        //! virtual copy constructor
        //!(needed for copying pointers to derived classes)
        virtual PointSurfaceSphere *Clone() const;


        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

//    protected:
        virtual const boost::scoped_ptr< DirectedPointSurface> &InitializeDefaultSphere( const float &DELTA);

    public:
        virtual std::vector< math::Vector3N> GetSurfacePointPositions() const;


        /////////////////////////
        //      READ/WRITE     //
        /////////////////////////


//        virtual std::ostream& Write( std::ostream &STREAM) const;
//
//        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM) const;
//
//        virtual std::istream& Read( std::istream &STREAM);

        virtual std::string GetClassName() const;

    }; // end class


    std::ostream &operator << ( std::ostream &STREAM, const PointSurfaceSphere &SPHERE);

    std::istream &operator >> ( std::istream &STREAM, PointSurfaceSphere &SPHERE);

} // end namespace geom




#endif /* POINT_SURFACE_SPHERE_H_ */

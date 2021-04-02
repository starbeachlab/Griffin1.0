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


#ifndef DIRECTED_SURFACE_POINT_H
#define DIRECTED_SURFACE_POINT_H

#include <deque>
#include <set>
#include <vector>
#include <cassert>
#include <iostream>

#include <boost/shared_ptr.hpp>

#include "surface_point.h"

namespace geom
{

    //! SurfacePoint contains information about an point representing the surface
    class DirectedSurfacePoint
    : public SurfacePoint
    {
     private:
         math::Vector3N  m_NormalVector;


     public:
         DirectedSurfacePoint()
         : SurfacePoint(),
         m_NormalVector()
         {}

        DirectedSurfacePoint( const DirectedSurfacePoint &SURF_POINT)
          :  SurfacePoint( SURF_POINT),
          m_NormalVector( SURF_POINT.m_NormalVector)
        {}

        DirectedSurfacePoint(  const float &X,  const float &Y, const float &Z,  const float &NX,  const float &NY, const float &NZ)
        : SurfacePoint( X, Y, Z),
        m_NormalVector( NX, NY, NZ)
        {}

        DirectedSurfacePoint( const math::Vector3N &POS, const math::Vector3N &NORM, const float &SURF = 0.0)
        : SurfacePoint( POS, SURF),
        m_NormalVector( NORM)
        { DebugWrite( __FUNCTION__);}

        virtual ~DirectedSurfacePoint(){}

        virtual DirectedSurfacePoint *Clone() const{ /*std::cout << "Directed Clone" << std::endl;*/ return new DirectedSurfacePoint( *this);}

        virtual void SetNormalVector( const float &X,  const float &Y, const float &Z){ m_NormalVector(0) = X; m_NormalVector(1) = Y; m_NormalVector(2) = Z;}
        virtual void SetNormalVector( const math::Vector3N &NORM){ m_NormalVector = NORM;}

        virtual math::Vector3N GetNormalVector() const{ return m_NormalVector;}



        virtual std::ostream& Write( std::ostream &STREAM) const
        {
//            std::cout << "directed WriteVmdCommands" << std::endl;

//            STREAM << GetClassName() << std::endl;
            STREAM << mystr::GetClassName( __PRETTY_FUNCTION__) << std::endl;
            SurfacePoint::Write( STREAM);
            STREAM << m_NormalVector << std::endl;
            return STREAM;
        }

        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM, const size_t &MOL_ID = 0) const
        {
            SurfacePoint::WriteVmdCommands( STREAM, MOL_ID);
//            STREAM << "graphics 0 color white" <<  std::endl;
            math::Vector3N end( SurfacePoint::GetPosition() + m_NormalVector);
//            math::Vector3N middle( SurfacePoint::GetPosition() + 0.9 * m_NormalVector);
//            STREAM << "graphics 0 cylinder {" << SurfacePoint::GetPosition()( 0) << " " << SurfacePoint::GetPosition()( 1) << " " << SurfacePoint::GetPosition()( 2) << "} {" << middle( 0) << " "  << middle( 1) << " "  << middle( 2) << "} radius 0.04" << std::endl;
//            STREAM << "graphics 0 cone {" << middle( 0) << " "  << middle( 1) << " "  << middle( 2) << "} {" << end( 0) << " "  << end( 1) << " "  << end( 2) << "} radius 0.06" << std::endl;
            STREAM << "graphics " << MOL_ID << " cone {" << SurfacePoint::GetPosition()( 0) << " " << SurfacePoint::GetPosition()( 1) << " " << SurfacePoint::GetPosition()( 2) << "} {" << end( 0) << " "  << end( 1) << " "  << end( 2) << "} radius 0.06" << std::endl;
            return STREAM;
        }

        virtual std::ostream &WriteAsPdb( std::ostream &STREAM, const size_t &MOL_ID = 0, const size_t &ATOM_ID = 0, const std::string ATOM = "H")
        {
          STREAM << "ATOM  ";
          STREAM.width( 5);
          STREAM <<  ATOM_ID % int( 1e5) << "  ";
          STREAM.width( 3);
          STREAM << ATOM << " SUR ";
          STREAM.width(5);
          STREAM << MOL_ID % int( 1e5) << "    ";
          STREAM.width( 8);
          STREAM.precision( 3);
          STREAM << SurfacePoint::GetPosition()( 0);
          STREAM.width( 8);
          STREAM.precision( 3);
          STREAM << SurfacePoint::GetPosition()( 1);
          STREAM.width( 8);
          STREAM.precision( 3);
          STREAM << SurfacePoint::GetPosition()( 2);
          STREAM.width( 15);
          STREAM << " ";
          STREAM << std::endl;
          return STREAM;
        }

        virtual std::istream &Read( std::istream &STREAM)
        {
            std::string str;
            STREAM >> str;
            assert( str == GetClassName());
            SurfacePoint::Read( STREAM);
            STREAM >> m_NormalVector;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
//            std::cout << __FUNCTION__ << std::endl;
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class DirectedSurfacePoint
} // end namespace geom


#endif

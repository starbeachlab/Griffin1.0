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


#ifndef CYLINDER_H_
#define CYLINDER_H_

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"

#include "object1D.h"

namespace geom
{

    class Cylinder
    : public Object1D
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        float m_Radius;
        float m_Height;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        Cylinder()
        : Object1D(),
        m_Radius(),
        m_Height()
        {}

        //! construct from position, radius and height
        Cylinder( const math::Vector3N &POSITION, const float &RADIUS, const float &HEIGHT)
        : Object1D( POSITION),
        m_Radius( RADIUS),
        m_Height( HEIGHT)
        {}

        //! construct from position, z axis, radius and height
        Cylinder( const math::Vector3N &POSITION, const math::Vector3N &Z_AXIS, const float &RADIUS, const float &HEIGHT)
        : Object1D( POSITION, Z_AXIS),
        m_Radius( RADIUS),
        m_Height( HEIGHT)
        {}

        //! copy constructor
        Cylinder( const Cylinder &ORIGINAL)
        : Object1D( ORIGINAL),
        m_Radius( ORIGINAL.m_Radius),
        m_Height( ORIGINAL.m_Height)
        {}

        //! virtual destructor
        virtual ~Cylinder(){}

        //! virtual copy constructor
        virtual Cylinder *Clone() const{ return new Cylinder( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual std::vector< float> GetSizes() const
        {
            std::vector< float> sizes( 2);
            sizes[ 0] = m_Radius;
            sizes[ 1] = m_Height;
            return sizes;
        }

        virtual std::vector< math::Vector3N> GetAxes() const
        {
            return Object1D::GetAxes();
        }

        virtual float GetTotalSurface() const
        {
            return 2.0 * math::Pi * m_Radius * ( m_Radius + m_Height);
        }

        virtual float MaxRadius() const
        { return std::sqrt( math::Square( m_Radius) + math::Square( 0.5 * m_Height));}

        virtual bool IsPointWithin( const math::Vector3N &POS) const
        {
            math::Vector3N
                connect( POS - GetPosition());
            float
                angle( math::Angle( connect, m_ZAxis)),
                dist( connect.Length());

            if
            (
                    dist * cos( angle) < m_Height
                    &&
                    dist * sin( angle) < m_Radius
            )
            {
                return true;
            }
            return false;
        }

        virtual store::Limits3D BoxLimits() const
        {
            math::Vector3N
                axis( 0.5 * m_Height * m_ZAxis);

            // these are estimations of the limits
            float
                xdist( math::ProjectionOnVector( axis, math::Vector3N( 1.0, 0.0, 0.0)) + m_Radius),
                ydist( math::ProjectionOnVector( axis, math::Vector3N( 0.0, 1.0, 0.0)) + m_Radius),
                zdist( math::ProjectionOnVector( axis, math::Vector3N( 0.0, 0.0, 1.0)) + m_Radius);
            return store::Limits3D
            (
                    m_Position[0] - xdist,
                    m_Position[0] + xdist,
                    m_Position[1] - ydist,
                    m_Position[1] + ydist,
                    m_Position[2] - zdist,
                    m_Position[2] + zdist
            );
        }



        // point on the surface that has shortest distance
        virtual math::Vector3N ProjectionOnSurface( const math::Vector3N &POS) const
        {
            math::Vector3N
                connect( POS - GetPosition()),
                xaxis( math::CrossProduct( math::CrossProduct( Object1D::m_ZAxis, connect), Object1D::m_ZAxis));

            float
                angle( math::Angle( m_ZAxis, connect)),
                length( connect.Length()),
                z( length * cos( angle));

            if( z < m_Height)
            {
                return m_Radius / tan( angle) * Object1D::m_ZAxis + m_Radius * xaxis;
            }

            float
                r( length * sin( angle));

            if( r < m_Radius)
            {
                return m_Height * Object1D::m_ZAxis + m_Height / tan( angle) * xaxis;
            }

            return m_Height * Object1D::m_ZAxis + m_Radius * xaxis;
        }



        virtual float ClosestDistance( const math::Vector3N &POS) const // TODO: QUATSCH!!!
        {
            math::Vector3N
                connect( POS - GetPosition());
            float
                angle( math::Angle( m_ZAxis, connect)),
                length( connect.Length()),
                z( length * cos( angle));

            if( z < m_Height)
            {
                return length * sin( angle) - m_Radius;
            }

            float
                r( length * sin( angle));

            if( r < m_Radius)
            {
                return length * cos( angle) - m_Height;
            }

            float
                diagonal( sqrt( math::Square( m_Height) + math::Square( m_Radius)));

            return math::LawOfCosinus( length, diagonal, angle - asin( m_Radius / diagonal));
        }


        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            return STREAM;
        }

        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM, const size_t &MOL_ID = 0) const
        {
            math::Vector3N
                pos1( GetPosition() - 0.5 * m_Height * m_ZAxis),
                pos2( GetPosition() + 0.5 * m_Height * m_ZAxis);
            STREAM << "graphics " << MOL_ID << " cylinder {" << pos1( 0) << " " << pos1( 1) << " " << pos1( 2) << "} {" << pos2( 0) << " " << pos2( 1) << " " << pos2( 2) << "} radius " << m_Radius << " resolution 30 filled yes" << std::endl;
            return STREAM;
        }

//        virtual std::string GetClassName() const
//        {
//            return mystr::GetClassName( __PRETTY_FUNCTION__);
//        }

    }; // end class Cylinder
} // end namespace geom




#endif /* CYLINDER_H_ */

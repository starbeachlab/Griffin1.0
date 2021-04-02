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


#ifndef GEOM_SPHERE_H
#define GEOM_SPHERE_H

#include "object.h"
#include "../math/vector3N.h"
#include "../math/constants.h"

namespace geom
{

    class Sphere
    : public Object
    {
    protected:
        float     m_Radius;
    public:
        //! The default constructor, initializes all values as zero.
        Sphere( const float RADIUS = 0.0)
        : Object( math::Vector3N( 0.0, 0.0, 0.0)),
        m_Radius( RADIUS)
        { /*std::cout << "Sphere: construct from radius: " << m_Radius << std::endl;*/}

        //! Single value constructor.
        Sphere( const float& X, const float& Y, const float& Z, const float& RADIUS)
        : Object( math::Vector3N( X, Y, Z)),
        m_Radius( RADIUS)
        {}

        //! Construct from position and radius
        Sphere( const math::Vector3N& POSITION, const float& RADIUS)
        : Object( POSITION),
        m_Radius( RADIUS)
        {
//            std::cout << "Sphere: construct from position and radius" << std::endl;
            assert( RADIUS > 0);
        }

        //! Copy constructor
        Sphere( const Sphere& SPHERE)
        : Object( SPHERE),
        m_Radius( SPHERE.m_Radius)
        { /*std::cout << "Sphere: copy constructor" << std::endl;*/}


        //! Virtual copy constructor
        virtual Sphere* Clone() const
        { return new Sphere( *this);}

        //! Virtual destructor.
        virtual ~Sphere(){}

        virtual float GetRadius() const
        { return m_Radius;}

        virtual float MaxRadius() const
        { return m_Radius;}

        virtual void SetRadius( const float& RADIUS)
        { m_Radius = RADIUS;}

        virtual store::Limits3D BoxLimits() const
        {
            return store::Limits3D
            (
                    m_Position[0] - m_Radius,
                    m_Position[0] + m_Radius,
                    m_Position[1] - m_Radius,
                    m_Position[1] + m_Radius,
                    m_Position[2] - m_Radius,
                    m_Position[2] + m_Radius
            );
        }



        virtual float TotalSurface() const
        {
            return 2.0 * math::Pi * math::Square( m_Radius);
        }

        virtual float Volume() const
        {
            return ( 4.0 / 3.0) * math::Pi * pow( m_Radius, 3);
        }

//        virtual math::Vector3N GetPosition() const{ return math::Vector3N( *this);}

//        virtual void SetPosition( const math::Vector3N &NEWPOS){ *this = NEWPOS;}

        virtual bool IsPointWithin( const math::Vector3N & POS) const
        {
            if( ( POS - GetPosition()).Length() < m_Radius)
            { return true;}
            return false;
        }

        virtual float ClosestDistance( const math::Vector3N &POS) const // TODO: QUATSCH!!!
        {
            return ( POS - GetPosition()).Length() - m_Radius;
        }

        virtual math::Vector3N ProjectionOnSurface( const math::Vector3N &POS) const
        {
            return GetPosition() + ( POS - GetPosition()).SetToLength( m_Radius);
        }

        virtual bool IsPointOnSurface( const math::Vector3N) const { return false;}

        virtual std::vector< float> GetSizes() const{ std::cout << "Sphere::GetSizes()" << std::endl; std::vector< float> v(1); v[ 0]= m_Radius; return v;}

        virtual std::ostream& Write( std::ostream &STREAM) const
        {
            STREAM << mystr::GetClassName( __PRETTY_FUNCTION__) << std::endl;
            STREAM << GetPosition()( 0)<< " " << GetPosition()( 1) << " " << GetPosition()( 2) << std::endl;
            STREAM << m_Radius << std::endl;
            return STREAM;
        }

        virtual std::istream& Read( std::istream &STREAM)
        {
//            STREAM << m_Position[0] << "  " << m_Position[1] << "  " << m_Position[2]<< "  " << m_NormalVector[0] << "  " << m_NormalVector[1] << "  " << m_NormalVector[2] << std::endl;
            return STREAM;
        }

        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM, const size_t &MOL_ID = 0) const
        {
    //        STREAM << "mol new" << std::endl;
            STREAM << "graphics " << MOL_ID << " color red" <<  std::endl;
            STREAM << "graphics " << MOL_ID << " sphere {" << GetPosition()(0) << " " << GetPosition()(1) << " " << GetPosition()(2) << "} radius " << m_Radius << " resolution 21" << std::endl;
            STREAM << "graphics " << MOL_ID << " sphere {" << GetPosition()(0) << " " << GetPosition()(1) << " " << GetPosition()(2) << "} radius 0.1 resolution 21" << std::endl;
            return STREAM;
        }

    }; // end class Sphere


    //! Generate random point on the sphere surface (used in MAYER)
    inline
    math::Vector3N GenerateRandomPointOnUnitSphereSurface()
    {
//        StandardWrite( __FUNCTION__);
        math::Vector3N point;
        float squared_length;
        do
        {
            point.Randomize( -1.0, 1.0);
            squared_length = point.SquaredLength();
        } while( squared_length == 0 || squared_length > 1.0);
        return point.Normalize();
    }



    inline
    std::ostream &operator << ( std::ostream &STREAM, const Sphere &SPHERE)
    { return SPHERE.Write( STREAM);}


} // end namespace geom

#endif

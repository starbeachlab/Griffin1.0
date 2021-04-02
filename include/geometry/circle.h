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


#ifndef CIRCLE_H_
#define CIRCLE_H_

#include "plane.h"



namespace geom
{

    class Circle
    : public Plane
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        float m_Radius;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        Circle()
        : Plane(),
        m_Radius()
        {}

        //! construct from data
        Circle( const math::Vector3N &POSITION, const math::Vector3N &NORMAL_VECTOR, const float &RADIUS)
        : Plane( POSITION, NORMAL_VECTOR),
        m_Radius( RADIUS)
        {}

        //! copy constructor
        Circle( const Circle &ORIGINAL)
        : Plane( ORIGINAL),
        m_Radius( ORIGINAL.m_Radius)
        {}

        //! virtual destructor
        virtual ~Circle(){}

        //! virtual copy constructor
        virtual Circle *Clone() const{ return new Circle( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual float GetRadius() const
        {
            return m_Radius;
        }

        virtual void SetRadius( const float &RADIUS)
        {
            m_Radius = RADIUS;
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
            Plane::Write( STREAM);
            STREAM << "circle_radius: " << m_Radius << std::endl;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class Circle
} // end namespace geom




#endif /* CIRCLE_H_ */

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


#ifndef TORUS_H_
#define TORUS_H_

#include "circle.h"


namespace geom
{

    class Torus
    : public Circle
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        float m_RingRadius;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        Torus()
        : Circle(),
        m_RingRadius()
        {}

        //! construct from data
        Torus
        (
                const math::Vector3N &POSITION,
                const math::Vector3N &NORMAL_VECTOR,
                const float &CIRCLE_RADIUS,
                const float &RING_RADIUS
        )
        : Circle( POSITION, NORMAL_VECTOR, CIRCLE_RADIUS),
        m_RingRadius( RING_RADIUS)
        {}

        //! copy constructor
        Torus( const Torus &ORIGINAL)
        : Circle( ORIGINAL),
        m_RingRadius( ORIGINAL.m_RingRadius)
        {}

        //! virtual destructor
        virtual ~Torus(){}

        //! virtual copy constructor
        virtual Torus *Clone() const{ return new Torus( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        const float &GetTorusRadius() const
        {
            return m_RingRadius;
        }

        void SetTorusRadius( const float &TORUS_RADIUS)
        {
            m_RingRadius = TORUS_RADIUS;
        }

        float GetMaxRadius() const
        {
            return m_RingRadius + Circle::m_Radius;
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
            Circle::Write( STREAM);
            STREAM << "torus_radius: "<< m_RingRadius << std::endl;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class Torus
} // end namespace geom




#endif /* TORUS_H_ */

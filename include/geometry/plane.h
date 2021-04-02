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


#ifndef PLANE_H_
#define PLANE_H_

#include "object.h"


namespace geom
{

    class Plane
    : public Object
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        math::Vector3N m_NormalVector;
        math::Vector3N m_ArbitraryXAxis;
        math::Vector3N m_ArbitraryYAxis;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        Plane()
        : Object(),
        m_NormalVector()
        {}

        //! construct from math::Vector3N
        Plane( const math::Vector3N &POSITION, const math::Vector3N &NORMAL_VECTOR)
        : Object( POSITION),
        m_NormalVector( NORMAL_VECTOR),
        m_ArbitraryXAxis(),
        m_ArbitraryYAxis()
        {}

        //! copy constructor
        Plane( const Plane &ORIGINAL)
        : Object( ORIGINAL),
        m_NormalVector( ORIGINAL.m_NormalVector),
        m_ArbitraryXAxis( ORIGINAL.m_ArbitraryXAxis),
        m_ArbitraryYAxis( ORIGINAL.m_ArbitraryYAxis)
        {}

        //! virtual destructor
        virtual ~Plane(){}

        //! virtual copy constructor
        virtual Plane *Clone() const{ return new Plane( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        const math::Vector3N &GetNormalVector() const
        {
            return m_NormalVector;
        }

        void SetNormalVector( const math::Vector3N &NORMAL_VECTOR)
        {
            m_NormalVector = NORMAL_VECTOR;
        }

        const math::Vector3N &GetArbitraryXAxis() const
        {
            return m_ArbitraryXAxis;
        }

        void SetArbitraryXAxis( const math::Vector3N &AXIS)
        {
            m_ArbitraryXAxis = AXIS;
        }

        const math::Vector3N &GetArbitraryYAxis() const
        {
            return m_ArbitraryYAxis;
        }

        void SetArbitraryYAxis( const math::Vector3N &AXIS)
        {
            m_ArbitraryYAxis = AXIS;
        }


        void CalculateArbitraryXAndYAxis()
        {
            if( std::abs( m_NormalVector( 1)) > 0.001)  // sort of a fast angle check..
            {
                m_ArbitraryXAxis = math::CrossProduct( math::Vector3N( 1.0, 0.0, 0.0), m_NormalVector).NormalizedCopy();
            }
            else
            {
                m_ArbitraryXAxis = math::CrossProduct( math::Vector3N( 0.0, 1.0, 0.0), m_NormalVector).NormalizedCopy();
            }
            m_ArbitraryYAxis = math::CrossProduct( m_NormalVector, m_ArbitraryXAxis).NormalizedCopy();

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
            STREAM << "positition: " << GetPosition() << std::endl;
            STREAM << "normal_vector: " << m_NormalVector << std::endl;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class Plane
} // end namespace geom




#endif /* PLANE_H_ */

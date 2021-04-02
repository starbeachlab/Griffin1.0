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


#ifndef TORUS_SLICE_H_
#define TORUS_SLICE_H_

#include "torus.h"
#include "../math/constants.h"

namespace geom
{

    class TorusSlice
    : public Torus
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        float m_MinimumAngle;
        float m_MaximumAngle;
        float m_DeltaOuterCircle;
        float m_DeltaInnerRing;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        TorusSlice()
        : Torus(),
        m_MinimumAngle(),
        m_MaximumAngle(),
        m_DeltaOuterCircle(),
        m_DeltaInnerRing()
        {}

        //! construct from data
        TorusSlice
        (
                const math::Vector3N &POSITION,
                const math::Vector3N &NORMAL_VECTOR,
                const float &CIRCLE_RADIUS,
                const float &RING_RADIUS,
                const float &MIN_ANGLE,
                const float &MAX_ANGLE,
                const float &DELTA_OUTER_CIRCLE,
                const float &DELTA_INNER_RING
        )
        : Torus( POSITION, NORMAL_VECTOR, CIRCLE_RADIUS, RING_RADIUS),
        m_MinimumAngle( MIN_ANGLE),
        m_MaximumAngle( MAX_ANGLE),
        m_DeltaOuterCircle( DELTA_OUTER_CIRCLE),
        m_DeltaInnerRing( DELTA_INNER_RING)
        {}

        //! copy constructor
        TorusSlice( const TorusSlice &ORIGINAL)
        : Torus( ORIGINAL),
        m_MinimumAngle( ORIGINAL.m_MinimumAngle),
        m_MaximumAngle( ORIGINAL.m_MaximumAngle),
        m_DeltaOuterCircle( ORIGINAL.m_DeltaOuterCircle),
        m_DeltaInnerRing( ORIGINAL.m_DeltaInnerRing)
        {}

        //! virtual destructor
        virtual ~TorusSlice(){}

        //! virtual copy constructor
        virtual TorusSlice *Clone() const{ return new TorusSlice( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        const float &GetMinimumAngle() const
        {
            return m_MinimumAngle;
        }

        void SetMinimumAngle( const float &MIN_ANGLE)
        {
            m_MinimumAngle = MIN_ANGLE;
        }

        const float &GetMaximumAngle() const
        {
            return m_MaximumAngle;
        }

        void SetMaximumAngle( const float &MAX_ANGLE)
        {
            m_MaximumAngle = MAX_ANGLE;
        }

        const float &GetDeltaOuterCircle() const
        {
            return m_DeltaOuterCircle;
        }

        void SetDeltaOuterCircle( const float &DELTA)
        {
            m_DeltaOuterCircle = DELTA;
        }

        const float &GetDeltaInnerRing() const
        {
            return m_DeltaInnerRing;
        }

        void SetDeltaInnerRing( const float &DELTA)
        {
            m_DeltaInnerRing = DELTA;
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
            Torus::Write( STREAM);
            STREAM << "minimum_angle: " << m_MinimumAngle << std::endl;
            STREAM << "maximum_angle: " << m_MaximumAngle << std::endl;
            STREAM << "delta_circle: " << m_DeltaOuterCircle << std::endl;
            STREAM << "delta_ring: " << m_DeltaInnerRing << std::endl;
            return STREAM;
        }

        virtual std::ostream &WriteVmdSurface( std::ostream &STREAM) const
        {
            // call CalculateArbitraryXAndYAxis(); before this!!
            math::Vector3N circle_point, torus_point, axis, norm;
            STREAM << "mol new" << std::endl;
            norm = GetPosition() + GetNormalVector();
            STREAM << "graphics 0 color white" << std::endl;
            STREAM << "graphics 0 cone {" << GetPosition()( 0) << " " << GetPosition()( 1) << " " << GetPosition()(2) << "} {" << norm( 0) << " " << norm( 1) << " " << norm( 2) << "} radius 0.1 resolution 3" << std::endl;
            STREAM << "graphics 0 color red" << std::endl;
            float max( m_MaximumAngle);
            if( m_MinimumAngle > max)
            {
                max += 2.0 * math::Pi;
            }
            for( float alpha = 0.0; alpha < 2 * math::Pi; alpha += m_DeltaOuterCircle)
            {
                circle_point = GetPosition() + GetRadius() * ( cos( alpha) * m_ArbitraryXAxis + sin( alpha) * m_ArbitraryYAxis);
                axis = ( GetPosition() - circle_point).NormalizedCopy();
//                std::cout << axis << "\n" << circle_point << std::endl;
                for( float beta = m_MinimumAngle; beta <= max; beta += m_DeltaInnerRing)
                {
                    torus_point = circle_point + Torus::m_RingRadius * ( cos( beta) * axis + sin( beta) * GetNormalVector());
                    STREAM << "graphics 0 sphere {" << torus_point( 0) << " " << torus_point(1) << " " << torus_point( 2) << "} radius 0.05 resolution 1" << std::endl;
                }
            }
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class TorusSlice
} // end namespace geom




#endif /* TORUS_SLICE_H_ */

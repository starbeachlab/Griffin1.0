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


#ifndef LIMITS3D_H_
#define LIMITS3D_H_

#include <cmath>

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"

#include "../math/vector3N.h"


namespace store
{

    class Limits3D
    : public StreamOperator
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        float m_XMin;  // TODO: test vectors for min and max instead!
        float m_XMax;
        float m_YMin;
        float m_YMax;
        float m_ZMin;
        float m_ZMax;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        Limits3D()
        : m_XMin( std::numeric_limits< float>::max()),
        m_XMax( -std::numeric_limits< float>::max()),
        m_YMin( std::numeric_limits< float>::max()),
        m_YMax( -std::numeric_limits< float>::max()),
        m_ZMin( std::numeric_limits< float>::max()),
        m_ZMax( -std::numeric_limits< float>::max())
        {}

        //! construct from data
        Limits3D
        (
                const float &XMIN,
                const float &XMAX,
                const float &YMIN,
                const float &YMAX,
                const float &ZMIN,
                const float &ZMAX
        )
        : m_XMin( XMIN),
        m_XMax( XMAX),
        m_YMin( YMIN),
        m_YMax( YMAX),
        m_ZMin( ZMIN),
        m_ZMax( ZMAX)
        {}

        //! copy constructor
        Limits3D( const Limits3D &ORIGINAL)
        : m_XMin( ORIGINAL.m_XMin),
        m_XMax( ORIGINAL.m_XMax),
        m_YMin( ORIGINAL.m_YMin),
        m_YMax( ORIGINAL.m_YMax),
        m_ZMin( ORIGINAL.m_ZMin),
        m_ZMax( ORIGINAL.m_ZMax)
        {}

        //! virtual destructor
        virtual ~Limits3D(){}

        //! virtual copy constructor
        virtual Limits3D *Clone() const;

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual const float &GetXMin() const;

        virtual void SetXMin( const float &XMAX);

        virtual const float &GetXMax() const;

        virtual void SetXMax( const float &XMAX);

        virtual const float &GetYMin() const;

        virtual void SetYMin( const float &YMIN);

        virtual const float &GetYMax() const;

        virtual void SetYMax( const float &YMAX);

        virtual const float &GetZMin() const;

        virtual void SetZMin( const float &ZMIN);

        virtual const float &GetZMax() const;

        virtual void SetZMax( const float &ZMAX);

        virtual float XDelta() const;

        virtual float YDelta() const;

        virtual float ZDelta() const;

        virtual float GetMin( const size_t ID) const;

        virtual math::Vector3N GetMin() const;

        virtual float GetMax( const size_t ID) const;

        virtual math::Vector3N GetMax() const;

        virtual float GetDelta( const size_t ID) const;

        virtual void SetMin( const size_t ID, const float &VALUE);

        virtual void SetMax( const size_t ID, const float &VALUE);

        virtual void AdjustToVoxelSize( const float &VOXELSIZE);

        virtual
        bool
        AreLimitsInAgreementWithVoxelSize( const float &VOXELSIZE, const float &PRECISION) const;


        virtual Limits3D &Merge( const Limits3D &LIMITS);

        virtual bool IsWithin( const math::Vector3N &POSITION) const;


        virtual float Volume() const;

        virtual bool AreDefined() const;

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM);

        virtual std::ostream &Write( std::ostream &STREAM) const;

        virtual std::string GetClassName() const;

    }; // end class Limits3D
} // end namespace store




#endif /* LIMITS3D_H_ */

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


#ifndef SURF_FORCE_H_
#define SURF_FORCE_H_

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"



namespace phys
{

    class SurfForce
    : public Force
    {
    protected:
        ///////////////
        //    float     //
        ///////////////
        boost::shared_ptr< store::Map< std::string, geom::PointSurface> >      m_SurfMap;
        float                                                                 m_Max;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        SurfForce()
        : Force(),
        m_Max()
        {}

        //! construct from float
        SurfForce( const float&WEIGHT, const float &MAXIMUM, const boost::shared_ptr< store::Map< std::string, geom::PointSurface> > &SURF)
        : Force( WEIGHT),
        m_SurfMap( SURF),
        m_Max( MAXIMUM)
        {}

        //! copy constructor
        SurfForce( const SurfForce &ORIGINAL)
        : Force( ORIGINAL),
        m_Max( ORIGINAL.m_Max)
        {}

        //! virtual destructor
        virtual ~SurfForce(){}

        //! virtual copy constructor
        virtual SurfForce *Clone() const{ return new SurfForce( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        math::Vector3N operator()( const math::Vector3N &POS)
        {

        }

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

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

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class SurfForce
} // end namespace phys




#endif /* SURF_FORCE_H_ */

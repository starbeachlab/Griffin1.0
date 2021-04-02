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


#ifndef SURF_ATOM_H_
#define SURF_ATOM_H_


#include "atom.h"
#include "../geometry/point_surface_sphere.h"


namespace mol
{

    class SurfAtom
    : public Atom // floats position
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        boost::shared_ptr< geom::PointSurfaceSphere> m_Surf;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        SurfAtom()
        : Atom(),
        m_Surf()
        { std::cout << __FUNCTION__ << std::endl;}

        //! construct from position and type
        SurfAtom( const math::Vector3N &POS, const std::string &TYPE, const float &RESOLUTION = 10.0, const float &FACTOR = 1.0, const float &OFFSET = 0.3)
        : Atom( TYPE, POS)//,
//        m_Surf( boost::shared_ptr< geom::PointSurfaceSphere>( new geom::PointSurfaceSphere( POS, ( FACTOR * ElementDataMap()( TYPE)->GetVanDerWaalsRadius()) + OFFSET, RESOLUTION)))
        { std::cout << __PRETTY_FUNCTION__ << std::endl;}

        SurfAtom( const Atom &ATOM, const float &RESOLUTION = 10.0, const float &FACTOR = 1.0, const float &OFFSET = 0.3)
        : Atom( ATOM),
        m_Surf( boost::shared_ptr< geom::PointSurfaceSphere>( new geom::PointSurfaceSphere( ATOM.GetPosition(), ( FACTOR * ATOM.GetVanDerWaalsRadius()) + OFFSET, RESOLUTION)))
        {
            DebugWrite( __FUNCTION__ << " atom: " << ATOM.GetType() << " radius: " << ATOM.GetVanDerWaalsRadius() << " resolution: " << RESOLUTION << " factor: " << FACTOR << " offset: " << OFFSET);
        }


        //! copy constructor
        SurfAtom( const SurfAtom &ORIGINAL)
        : Atom( ORIGINAL),
        m_Surf( ORIGINAL.m_Surf)
        {}

        //! virtual destructor
        virtual ~SurfAtom(){}

        //! virtual copy constructor
        virtual SurfAtom *Clone() const{ return new SurfAtom( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual const boost::shared_ptr< geom::PointSurfaceSphere> &GetSurf() const
        { return m_Surf;}

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
            Atom::Read( STREAM);
            STREAM >> m_Surf;
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            Atom::Write( STREAM);
            STREAM << m_Surf;
            return STREAM;
        }

        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM) const
        {
//            Atom::WriteVmdCommands( STREAM);
            m_Surf->WriteVmdCommands( STREAM);
            return STREAM;
        }

    }; // end class Atom
} // end namespace mol




#endif /* SURF_ATOM_H_ */

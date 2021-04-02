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


#include "../../include/phys/potential_force.h"

namespace phys
{

    //! virtual copy constructor
    PotentialForce *PotentialForce::Clone() const{ return new PotentialForce( *this);}

    /////////////////////////
    //     OPERATORS       //
    /////////////////////////

    //! returns the potential force without requiring any actual information from atom; this function needs to be overwritten in derived classes for interactions requiring specific information of the atom
    math::Vector3N PotentialForce::operator()( const mol::Atom &ATOM) const
    { return m_Force;}


    float PotentialForce::Energy( const mol::Atom &ATOM) const
    {
    	return m_Energy;
    }

    void PotentialForce::operator = ( const math::Vector3N &FORCE)
    { m_Force = FORCE;}


//    void PotentialForce::operator += ( const math::Vector3N &FORCE, const float &ENERGY)
//    {
//        m_Force += FORCE;
//        m_ENERGY += ENERGY;
//    }


    /////////////////////////
    //      FUNCTIONS      //
    /////////////////////////

    void PotentialForce::AddToPotentialForce( const math::Vector3N &FORCE, const float &ENERGY)
    {
        m_Force += FORCE;
        m_Energy += ENERGY;
    }

    void PotentialForce::AddToPotentialForce( const math::Vector3N &FORCE)
    {
        m_Force += FORCE;
    }

    void PotentialForce::AddToPotentialForce( const float &ENERGY)
    {
         m_Energy += ENERGY;
    }

    /////////////////////////
    //      Read/Write     //
    /////////////////////////

    std::istream &PotentialForce::Read( std::istream &STREAM)
    {
       	DebugWrite( __PRETTY_FUNCTION__);
    	std::string str;
        STREAM >> m_Force[0] >> m_Force[1] >> m_Force[2];
        STREAM >> str;
        STREAM >> m_Energy;
        return STREAM;
    }

    std::ostream &PotentialForce::Write( std::ostream &STREAM) const
    {
        STREAM << GetClassName() << std::endl;
        STREAM << m_Force;
        STREAM << "energy: " << m_Energy << std::endl;
        return STREAM;
    }

    //! writes forces without class identifier; to be called from derived classes
    std::ostream &PotentialForce::WriteData( std::ostream &STREAM) const
    {
        STREAM << m_Force( 0) << " " << m_Force( 1) << "  " << m_Force( 2) << std::endl;
        STREAM << "energy: " << m_Energy << std::endl;
        return STREAM;
    }

    std::string PotentialForce::GetClassName() const
    {
        return mystr::GetClassName( __PRETTY_FUNCTION__);
    }

} // end namespace $

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


#ifndef ELECTROSTATIC_POTENTIAL_FORCE_H_
#define ELECTROSTATIC_POTENTIAL_FORCE_H_


#include "potential_force.h"
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"

namespace phys
{

    class PotentialElectrostaticForce
    : public PotentialForce
    {

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor, zero force
        PotentialElectrostaticForce()
        : PotentialForce()
        { /*std::cout << __PRETTY_FUNCTION__ << std::endl;*/}

        //! construct from force
        PotentialElectrostaticForce( const math::Vector3N &FORCE, const float &ENERGY)
        : PotentialForce( FORCE, ENERGY)
        {}

        //! copy constructor
        PotentialElectrostaticForce( const PotentialElectrostaticForce &ORIGINAL)
        : PotentialForce( ORIGINAL)
        {}

        //! virtual destructor
        virtual ~PotentialElectrostaticForce()
        {}

        //! virtual copy constructor
        virtual PotentialElectrostaticForce *Clone() const
        { return new PotentialElectrostaticForce( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        virtual math::Vector3N operator() ( const mol::Atom &ATOM) const
        {
            DebugWrite( __PRETTY_FUNCTION__);
            return ATOM.GetPartialCharge() * PotentialForce::m_Force;
        }

        virtual float Energy( const mol::Atom &ATOM) const
        {
            DebugWrite( __PRETTY_FUNCTION__);
            DebugWrite(  "coulomb: " << ATOM.GetPartialCharge() * PotentialForce::m_Energy << " charge: " << ATOM.GetPartialCharge());
        	return ATOM.GetPartialCharge() * PotentialForce::m_Energy;
        }


        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
           	DebugWrite( __PRETTY_FUNCTION__);
        	std::string str;
        	STREAM >> str;
        	if( str != GetClassName())
        	{
        		std::cout << "===> not the correct id <" << str << "> in " << __PRETTY_FUNCTION__ << std::endl;
        		exit( -1);
        	}
            PotentialForce::Read( STREAM);
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            STREAM << GetClassName() << "  ";
            PotentialForce::WriteData( STREAM);
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class PotentialElectrostaticForce
} // end namespace phys




#endif /* ELECTROSTATIC_POTENTIAL_FORCE_H_ */

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


#include "../../include/phys/potential_repulsive_vdw_force.h"



namespace phys
{
	math::Vector3N PotentialRepulsiveVanDerWaalsForce::s_Local = math::Vector3N();

	//! default constructor
	PotentialRepulsiveVanDerWaalsForce::PotentialRepulsiveVanDerWaalsForce()
	: PotentialForce()
	{}

	//! construct from data
	PotentialRepulsiveVanDerWaalsForce::PotentialRepulsiveVanDerWaalsForce( const math::Vector3N &FORCE, const float &ENERGY)
	: PotentialForce( FORCE, ENERGY)
	{}

	//! copy constructor
	PotentialRepulsiveVanDerWaalsForce::PotentialRepulsiveVanDerWaalsForce( const PotentialRepulsiveVanDerWaalsForce &ORIGINAL)
	: PotentialForce( ORIGINAL)
	{}

	//! virtual destructor
	PotentialRepulsiveVanDerWaalsForce::~PotentialRepulsiveVanDerWaalsForce(){}

	//! virtual copy constructor
	PotentialRepulsiveVanDerWaalsForce *PotentialRepulsiveVanDerWaalsForce::Clone() const{ return new PotentialRepulsiveVanDerWaalsForce( *this);}

	/////////////////////////
	//     OPERATORS       //
	/////////////////////////

	math::Vector3N PotentialRepulsiveVanDerWaalsForce::operator() ( const mol::Atom &ATOM) const
	{

#ifdef DEBUG
		std::cout << __PRETTY_FUNCTION__ << std::endl;
		std::cout << "epsilon: " << ATOM.GetEpsilon() << std::endl;
		std::cout << "vdw-radius: " << ATOM.GetVanDerWaalsRadius() << std::endl;
		std::cout << "potential_repulsive_vdw_force: " << ( PotentialForce::m_Force)( 0) << " " << ( PotentialForce::m_Force)( 1) << " " << ( PotentialForce::m_Force)( 2) << std::endl;
#endif


		s_Local = ATOM.GetRepulsiveVdwFactor() * PotentialForce::m_Force;
		return s_Local;
	}


	float PotentialRepulsiveVanDerWaalsForce::Energy( const mol::Atom &ATOM) const
	{
		DebugWrite( __PRETTY_FUNCTION__);
		DebugWrite(  "vdw: " << ATOM.GetRepulsiveVdwFactor()  * PotentialForce::m_Energy
				<< " epsilon: " << ATOM.GetEpsilon() << " radius: " << ATOM.GetVanDerWaalsRadius());

		return ATOM.GetRepulsiveVdwFactor() * PotentialForce::m_Energy;
	}
        /////////////////////////
        //      Read/Write     //
        /////////////////////////

	std::istream &PotentialRepulsiveVanDerWaalsForce::Read( std::istream &STREAM)
	{
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

	std::ostream &PotentialRepulsiveVanDerWaalsForce::Write( std::ostream &STREAM) const
	{
		STREAM << GetClassName() << " ";
		PotentialForce::WriteData( STREAM);
		return STREAM;
	}

	std::string PotentialRepulsiveVanDerWaalsForce::GetClassName() const
	{
		return mystr::GetClassName( __PRETTY_FUNCTION__);
	}




} // end namespace phys




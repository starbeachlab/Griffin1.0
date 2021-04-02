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


#include "../../include/phys/potential_attractive_vdw_force.h"



namespace phys
{
	math::Vector3N PotentialAttractiveVanDerWaalsForce::s_Local = math::Vector3N();

	//! default constructor
	PotentialAttractiveVanDerWaalsForce::PotentialAttractiveVanDerWaalsForce()
	: PotentialForce()
	{}

	//! construct from data
	PotentialAttractiveVanDerWaalsForce::PotentialAttractiveVanDerWaalsForce( const math::Vector3N &FORCE, const float &ENERGY)
	: PotentialForce( FORCE, ENERGY)
	{}

	//! copy constructor
	PotentialAttractiveVanDerWaalsForce::PotentialAttractiveVanDerWaalsForce( const PotentialAttractiveVanDerWaalsForce &ORIGINAL)
	: PotentialForce( ORIGINAL)
	{}

	//! virtual destructor
	PotentialAttractiveVanDerWaalsForce::~PotentialAttractiveVanDerWaalsForce(){}

	//! virtual copy constructor
	PotentialAttractiveVanDerWaalsForce *PotentialAttractiveVanDerWaalsForce::Clone() const{ return new PotentialAttractiveVanDerWaalsForce( *this);}

	/////////////////////////
	//     OPERATORS       //
	/////////////////////////

	math::Vector3N PotentialAttractiveVanDerWaalsForce::operator() ( const mol::Atom &ATOM) const
	{

#ifdef DEBUG
		std::cout << __PRETTY_FUNCTION__ << std::endl;
		std::cout << "epsilon: " << ATOM.GetEpsilon() << std::endl;
		std::cout << "vdw-radius: " << ATOM.GetVanDerWaalsRadius() << std::endl;
		std::cout << "potential_attractive_vdw_force: " << ( PotentialForce::m_Force)( 0) << " " << ( PotentialForce::m_Force)( 1) << " " << ( PotentialForce::m_Force)( 2) << std::endl;
#endif

//            return
//                math::Vector3N
//                (
//                        new math::Vector3N( sqrt( fabs( ATOM.GetEpsilon())) * pow( ATOM.GetVanDerWaalsRadius(), 3) * *PotentialForce::m_Force)  // TODO !!
//                );

		return s_Local = ATOM.GetAttractiveVdwFactor() * PotentialForce::m_Force;
	}


	float PotentialAttractiveVanDerWaalsForce::Energy( const mol::Atom &ATOM) const
	{
		DebugWrite( __PRETTY_FUNCTION__);
		DebugWrite(  "vdw: " << ATOM.GetAttractiveVdwFactor()  * PotentialForce::m_Energy
				<< " epsilon: " << ATOM.GetEpsilon() << " radius: " << ATOM.GetVanDerWaalsRadius());

//        	return sqrt( fabs( ATOM.GetEpsilon())) * pow( ATOM.GetVanDerWaalsRadius(), 3) * PotentialForce::m_Energy;

		return ATOM.GetAttractiveVdwFactor() * PotentialForce::m_Energy;
	}
        /////////////////////////
        //      Read/Write     //
        /////////////////////////

	std::istream &PotentialAttractiveVanDerWaalsForce::Read( std::istream &STREAM)
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

	std::ostream &PotentialAttractiveVanDerWaalsForce::Write( std::ostream &STREAM) const
	{
		STREAM << GetClassName() << " ";
		PotentialForce::WriteData( STREAM);
		return STREAM;
	}

	std::string PotentialAttractiveVanDerWaalsForce::GetClassName() const
	{
		return mystr::GetClassName( __PRETTY_FUNCTION__);
	}




} // end namespace phys




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


#include "../../include/phys/electrostatic_force.h"


namespace phys
{
	float ElectrostaticForce::sm_ForceEnergyUnitFactor = 332.0716;

    //! default constructor
    ElectrostaticForce::ElectrostaticForce( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE, const float &CUTOFF)
    : Force( CUTOFF),
    m_Molecule( IMPLICIT_MOLECULE)
    {
    	StandardWrite( __FUNCTION__ << " cutoff: " << CUTOFF);
    }

    ElectrostaticForce::ElectrostaticForce( const float &CUTOFF)
    : Force( CUTOFF),
    m_Molecule()
    { DebugWrite( __PRETTY_FUNCTION__);}


    //! copy constructor
    ElectrostaticForce::ElectrostaticForce( const ElectrostaticForce &ORIGINAL)
    : Force( ORIGINAL),
    m_Molecule( ORIGINAL.m_Molecule)
    {}

    //! virtual destructor
    ElectrostaticForce::~ElectrostaticForce(){}

    //! virtual copy constructor
    ElectrostaticForce *ElectrostaticForce::Clone() const{ return new ElectrostaticForce( *this);}

    /////////////////////////
    //     OPERATORS       //
    /////////////////////////

    //! the force caused by ATOM_B acting on ATOM_A
    math::Vector3N ElectrostaticForce::operator()( const mol::Atom &ATOM_A, const mol::Atom &ATOM_B) const
    {
        math::Vector3N connector( ATOM_A.GetPosition() - ATOM_B.GetPosition());
        float distance( connector.Length());
        if( distance > m_Cutoff)
        { return math::Vector3N();}
        return
			sm_ForceEnergyUnitFactor
			* ATOM_A.GetPartialCharge()
			* ATOM_B.GetPartialCharge()
			* std::pow( distance, -3)
			* connector;  // dist**-3 because connector is not normalized!
    }

    /////////////////////////
    //      FUNCTIONS      //
    /////////////////////////

    boost::shared_ptr< phys::PotentialForce> ElectrostaticForce::PotentialForce( const math::Vector3N & POS) const
    {
         DebugWrite( __PRETTY_FUNCTION__ );
        math::Vector3N
            potential_force;
        math::Vector3N
            connector;
        float
            distance,
            energy = 0.0;

        for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator mol_itr = m_Molecule->GetAtoms().begin(); mol_itr != m_Molecule->GetAtoms().end(); ++mol_itr)
        {
            connector = POS - ( *mol_itr)->GetPosition();
            distance = connector.Length();

#ifdef DEBUG
//            std::cout << "type: " << ( *mol_itr)->GetType() << " partial_charge: " << ( *mol_itr)->GetPartialCharge() << " distance: " << distance << " std::pow( distance, -3) * 332.0716 * partial_charge: <" << 332.0716 * ( *mol_itr)->GetPartialCharge() * std::pow( distance, -3) << "> (" << std::pow( distance, -3) << ") cutoff: " << m_Cutoff << std::endl;
#endif

            potential_force += ( *mol_itr)->GetPartialCharge() * std::pow( distance, -3) * connector;   // dist**-3 because connector is not normalized!
            energy += ( *mol_itr)->GetPartialCharge() / distance;

        }
        potential_force *= sm_ForceEnergyUnitFactor;
        energy *= sm_ForceEnergyUnitFactor;

#ifdef DEBUG
        std::cout << "potential_force: " << potential_force.Length()  << std::endl;
#endif

        return boost::shared_ptr< phys::PotentialForce>( new PotentialElectrostaticForce( potential_force, energy));
    }

    boost::shared_ptr< phys::PotentialForce> ElectrostaticForce::GetAssociatedPotentialForceObject() const
    {
    	/*std::cout << __PRETTY_FUNCTION__ << std::endl;*/
    	return boost::shared_ptr< phys::PotentialForce>( new PotentialElectrostaticForce());
    }


    void ElectrostaticForce::SetKjouleNmUnits()
    {
//    	sm_ForceEnergyUnitFactor = 138.935485; // 332.0716 * 0.4184 // Gromacs manual: unit:   kjoule / mol * nm * e^2
    	sm_ForceEnergyUnitFactor = 13893.5485; // 332.0716 * 41.84 // Gromacs manual: unit:   kjoule / mol * nm * e^2
    }


    /////////////////////////
    //      Read/Write     //
    /////////////////////////


    std::istream &ElectrostaticForce::Read( std::istream &STREAM)
    {
        return STREAM;
    }

    std::ostream &ElectrostaticForce::Write( std::ostream &STREAM) const
    {
        STREAM << GetClassName() << std::endl;
        return STREAM;
    }

    std::string ElectrostaticForce::GetClassName() const
    {
        return mystr::GetClassName( __PRETTY_FUNCTION__);
    }

} // end namespace phys


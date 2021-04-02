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


#include "../../include/phys/repulsive_vdw_force.h"


namespace phys
{
//		float RepulsiveVanDerWaalsForce::sm_ForceUnitFactor = 49152.0; // -( (-12) * 2**12);
//		float RepulsiveVanDerWaalsForce::sm_EnergyUnitFactor = 4096.0;  // 2**12;


        //! default constructor
        RepulsiveVanDerWaalsForce::RepulsiveVanDerWaalsForce( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE, const float &CUTOFF)
        : Force( CUTOFF),
        m_Molecule( IMPLICIT_MOLECULE)
        {
        	StandardWrite( __FUNCTION__ << " cutoff: " << CUTOFF);
        }

        //! default constructor
        RepulsiveVanDerWaalsForce::RepulsiveVanDerWaalsForce(  const float &CUTOFF)
        : Force( CUTOFF),
        m_Molecule()
        { DebugWrite( __PRETTY_FUNCTION__);}

        //! copy constructor (meaningless)
        RepulsiveVanDerWaalsForce::RepulsiveVanDerWaalsForce( const RepulsiveVanDerWaalsForce &ORIGINAL)
        : Force( ORIGINAL),
        m_Molecule( ORIGINAL.m_Molecule)
        {}

        //! virtual destructor
        RepulsiveVanDerWaalsForce::~RepulsiveVanDerWaalsForce(){}

        //! virtual copy constructor
        RepulsiveVanDerWaalsForce *RepulsiveVanDerWaalsForce::Clone() const{ return new RepulsiveVanDerWaalsForce( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////
        math::Vector3N RepulsiveVanDerWaalsForce::operator()( const mol::Atom &ATOM_A, const mol::Atom &ATOM_B) const
        {
            math::Vector3N connector = ATOM_A.GetPosition() - ATOM_B.GetPosition();
            float distance( connector.Length());
            connector.Normalize();
            if( distance > m_Cutoff)
            { return math::Vector3N();}
            return
				49152.0    // -( (-12) * 2**12);
				* sqrt( ATOM_A.GetEpsilon() * ATOM_B.GetEpsilon())
				* std::pow( ATOM_A.GetVanDerWaalsRadius() * ATOM_B.GetVanDerWaalsRadius(), 6)
				* std::pow( distance, -13)
				* connector;
            // sigma = 2 * vdw-radius // 12 * (2*2)**3 = 768
        }



        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////


        boost::shared_ptr< phys::PotentialForce>  RepulsiveVanDerWaalsForce::PotentialForce( const math::Vector3N & POSITION) const
        {
            DebugWrite( __PRETTY_FUNCTION__);
            math::Vector3N
                potential_force;
            math::Vector3N
                connector;
            float
                distance,
                energy = 0.0;

            // here we iterate through the molecule and calculate the potential force at POSITION
            for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator mol_itr = m_Molecule->GetAtoms().begin(); mol_itr != m_Molecule->GetAtoms().end(); ++mol_itr)
            {
                connector = POSITION - ( *mol_itr)->GetPosition();
                distance = connector.Length();

#ifdef DEBUG
//                std::cout << "type: " << ( *mol_itr)->GetType()
//                    << " epsilon: " << ( *mol_itr)->GetEpsilon()
//                    << " radius(minimum of lennard-jones): " << ( *mol_itr)->GetVanDerWaalsRadius()
//                    << " distance: " << distance
//                    << " | 49152.0 * sqrt( Epsilon) * std::pow( VDWRadius, 6) * std::pow( distance, -13): <" << 49152.0 * sqrt( fabs( ( *mol_itr)->GetEpsilon())) * std::pow( ( *mol_itr)->GetVanDerWaalsRadius(), 6) * std::pow( distance, -13)
//                    << "> (" << std::pow( distance, -13) << ") cutoff: " << m_Cutoff << std::endl;
#endif

                potential_force += ( *mol_itr)->GetRepulsiveVdwFactor() * std::pow( distance, -14) * connector;  // dist**-14 because connector is not normalized
                energy += ( *mol_itr)->GetRepulsiveVdwFactor() * std::pow( distance, -12);
            }
            potential_force *= 49152.0; // -( (-12) * 2**12)
            energy *= 4096.0;  // 2**12;

#ifdef DEBUG
            std::cout << "potential_force: " << potential_force.Length() << std::endl;
#endif

            return boost::shared_ptr< phys::PotentialForce>( new PotentialRepulsiveVanDerWaalsForce( potential_force, energy));;
        }

        //math::Vector3N PotentialForce( const math::Vector3N & POS, const boost::shared_ptr< mol::MoleculeIF> &MOLECULE);


        boost::shared_ptr< phys::PotentialForce> RepulsiveVanDerWaalsForce::GetAssociatedPotentialForceObject() const
        { return boost::shared_ptr< phys::PotentialForce>( new PotentialRepulsiveVanDerWaalsForce());}


//        void RepulsiveVanDerWaalsForce::SetKjouleNmUnits()
//        {
////        	sm_ForceUnitFactor = ;
////        	sm_EnergyUnitFactor = ;
//        }


        /////////////////////////
        //      Read/Write     //
        /////////////////////////
        std::istream &RepulsiveVanDerWaalsForce::Read( std::istream &STREAM)
        {
//            Force::Read( STREAM);
            return STREAM;
        }

        std::ostream &RepulsiveVanDerWaalsForce::Write( std::ostream &STREAM) const
        {
            STREAM << GetClassName();
//            Force::Write( STREAM);
            return STREAM;
        }



        std::string RepulsiveVanDerWaalsForce::GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

} // end namespace phys


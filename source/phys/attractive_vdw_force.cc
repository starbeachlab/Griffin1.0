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


#include "../../include/phys/attractive_vdw_force.h"


namespace phys
{
//		float AttractiveVanDerWaalsForce::sm_ForceUnitFactor = -768.0;     // -(-2 * (-6) * 2**6)
//		float AttractiveVanDerWaalsForce::sm_EnergyUnitFactor = -128.0;    // -2 * 2**6


        //! default constructor
        AttractiveVanDerWaalsForce::AttractiveVanDerWaalsForce( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE, const float &CUTOFF)
        : Force( CUTOFF),
        m_Molecule( IMPLICIT_MOLECULE)
        {
        	StandardWrite( __FUNCTION__ << " cutoff: " << CUTOFF);
        }

        //! default constructor
        AttractiveVanDerWaalsForce::AttractiveVanDerWaalsForce(  const float &CUTOFF)
        : Force( CUTOFF),
        m_Molecule()
        { DebugWrite( __PRETTY_FUNCTION__);}

        //! copy constructor (meaningless)
        AttractiveVanDerWaalsForce::AttractiveVanDerWaalsForce( const AttractiveVanDerWaalsForce &ORIGINAL)
        : Force( ORIGINAL),
        m_Molecule( ORIGINAL.m_Molecule)
        {}

        //! virtual destructor
        AttractiveVanDerWaalsForce::~AttractiveVanDerWaalsForce(){}

        //! virtual copy constructor
        AttractiveVanDerWaalsForce *AttractiveVanDerWaalsForce::Clone() const{ return new AttractiveVanDerWaalsForce( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////
        math::Vector3N AttractiveVanDerWaalsForce::operator()( const mol::Atom &ATOM_A, const mol::Atom &ATOM_B) const
        {
            math::Vector3N connector( ATOM_A.GetPosition() - ATOM_B.GetPosition());
            float distance( connector.Length());
            connector.Normalize();
            if( distance > m_Cutoff)
            { return math::Vector3N();}
            return
            // sigma = 2 * vdw-radius // 12 * (2*2)**3 = 768
				-768.0	     // -(-2 * (-6) * 2**6)
				* sqrt( ATOM_A.GetEpsilon() * ATOM_B.GetEpsilon())
				* std::pow( ATOM_A.GetVanDerWaalsRadius() * ATOM_B.GetVanDerWaalsRadius(), 3)
				* std::pow( distance, -7)
				* connector;
        }



        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        boost::shared_ptr< phys::PotentialForce>  AttractiveVanDerWaalsForce::PotentialForce( const math::Vector3N & POSITION) const
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
                distance = connector.Length();  // TODO : use squared length and adjust pow below!

#ifdef DEBUG
//                std::cout << "type: " << ( *mol_itr)->GetType()
//                    << " epsilon: " << ( *mol_itr)->GetEpsilon()
//                    << " radius(minimum of lennard-jones): " << ( *mol_itr)->GetVanDerWaalsRadius()
//                    << " distance: " << distance
//                    << " | -768.0 * sqrt( Epsilon) * std::pow( VDWRadius, 3) * std::pow( distance, -7): <" << -768.0 * sqrt( fabs( ( *mol_itr)->GetEpsilon())) * std::pow( ( *mol_itr)->GetVanDerWaalsRadius(), 3) * std::pow( distance, -7)
//                    << "> (" << std::pow( distance, -7) << ") cutoff: " << m_Cutoff << std::endl;
#endif

                potential_force += ( *mol_itr)->GetAttractiveVdwFactor() * std::pow( distance, -8) * connector;
                energy += ( *mol_itr)->GetAttractiveVdwFactor() * std::pow( distance, -6);
            }

            potential_force *= -768.0;     // -(-2 * (-6) * 2**6)  // see namd manual for lennard jones potential
            energy *= -128.0;              // -2 * 2**6


#ifdef DEBUG
            std::cout << "potential_force: " << potential_force.Length() << std::endl;
#endif

            return boost::shared_ptr< phys::PotentialForce>( new PotentialAttractiveVanDerWaalsForce( potential_force, energy));;
        }

        //math::Vector3N PotentialForce( const math::Vector3N & POS, const boost::shared_ptr< mol::MoleculeIF> &MOLECULE);


        boost::shared_ptr< phys::PotentialForce> AttractiveVanDerWaalsForce::GetAssociatedPotentialForceObject() const
        { return boost::shared_ptr< phys::PotentialForce>( new PotentialAttractiveVanDerWaalsForce());}

//        void AttractiveVanDerWaalsForce::SetKjouleNmUnits()
//        {
////        	sm_ForceUnitFactor = ;
////        	sm_EnergyUnitFactor = ;
//        }



        /////////////////////////
        //      Read/Write     //
        /////////////////////////
        std::istream &AttractiveVanDerWaalsForce::Read( std::istream &STREAM)
        {
//            Force::Read( STREAM);
            return STREAM;
        }

        std::ostream &AttractiveVanDerWaalsForce::Write( std::ostream &STREAM) const
        {
            STREAM << GetClassName();
//            Force::Write( STREAM);
            return STREAM;
        }



        std::string AttractiveVanDerWaalsForce::GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

} // end namespace phys


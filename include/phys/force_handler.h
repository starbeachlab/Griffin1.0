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


#ifndef FORCE_HANDLER_H_
#define FORCE_HANDLER_H_

namespace phys
{
    class ForceHandler
    {

        ForceHandler
        (
                const boost::shared_ptr<math::Function< float, float> > &PAIR_POTENTIAL,
                const boost::shared_ptr<math::Function< float, float> > &BACKBONE_DISTANCE_POTENTIAL,
                const boost::shared_ptr<math::Function< float, float> > &BACKBONE_ANGLE_POTENTIAL = boost::shared_ptr< math::Function< float, float> >(),
                const boost::shared_ptr<math::Function< float, float> > &BACKBONE_DIHEDRAL_POTENTIAL = boost::shared_ptr< math::Function< float, float> >(),
                const size_t &OMIT_NEXT_NEIGHBORS = 1,
                const bool &PAIR_SWITCH = false,
                const bool &ANGLE_SWITCH = false
        )
        : m_PairPotential( PAIR_POTENTIAL),
        m_DistancePotential( BACKBONE_DISTANCE_POTENTIAL),
        m_AnglePotential( BACKBONE_ANGLE_POTENTIAL),
        m_DihedralPotential( BACKBONE_DIHEDRAL_POTENTIAL),
        m_OmitedNeighbors( OMIT_NEXT_NEIGHBORS),
        m_IsPairSwitch( SWITCH),
        m_IsAngleSwitch( ANGLE_SWITCH)
        {}



        void operator()( const boost::shared_ptr< mol::SimpleMolecule< mol::AtomUnit> > &MOLECULE)
        {
            // pair forces
            for( std::vector< boost::shared_ptr< mol::AtomUnit> >::const_iterator first_atom_itr = MOLECULE->GetAtoms().begin(); first_atom_itr != MOLELECULE->GetAtoms().end(); ++first_atom_itr)
            {

                for( std::vector< boost::shared_ptr< mol::AtomUnit> >::const_iterator second_atom_itr = MOLECULE->GetAtoms().begin(); second_atom_itr <= first_atom_itr + m_OmitedNeighbors; ++second_atom_itr)
                {
                    m_PairForces->operator()( *first_atom_itr, *second_atom_itr, m_PairPotential);
                }

                for( std::vector< boost::shared_ptr< mol::AtomUnit> >::const_iterator second_atom_itr = first_atom_itr + m_OmitedNeighbors; second_atom_itr != MOLECULE->GetAtoms().end(); ++second_atom_itr)
                {
                    m_PairForces->operator()( *first_atom_itr, *second_atom_itr, m_PairPotential);
                }
            }

            if( m_IsSwitch)
            {
                m_PairPotential->Increment();
            }

            // backbone distances
            for( std::vector< boost::shared_ptr< mol::AtomUnit> >::const_iterator atom_itr = MOLECULE->GetAtoms().begin(); atom_itr + 1 != MOLELECULE->GetAtoms().end(); ++atom_itr)
            {
                m_DistancePotential->operator()( atom_itr);
            }

            if( m_AnglePotential != NULL)
            {
                // backbone angles
                for( std::vector< boost::shared_ptr< mol::AtomUnit> >::const_iterator atom_itr = MOLECULE->GetAtoms().begin(); atom_itr + 2 != MOLELECULE->GetAtoms().end(); ++atom_itr)
                {
                    m_AnglePotential->operator()( atom_itr);
                }
                if( m_IsAngleSwitch)
                {
                    m_AnglePotential->Increment();
                }
            }
            if( m_Dihedral_Potential != NULL)
            {
                // backbone dihedrals
                for( std::vector< boost::shared_ptr< mol::AtomUnit> >::const_iterator atom_itr = MOLECULE->GetAtoms().begin(); atom_itr + 3 != MOLELECULE->GetAtoms().end(); ++atom_itr)
                {
                    m_DihedralPotential->operator()( atom_itr);
                }
                if( m_IsAngleSwitch)
                {
                    m_DihedralPotential->Increment();
                }
            }
        }
    };
}


#endif /* FORCE_HANDLER_H_ */

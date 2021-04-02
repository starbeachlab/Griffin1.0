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


#ifndef VDW_ENERGY_H_
#define VDW_ENERGY_H_

#include <>

#include "../readwrite/class_name.h"


namespace phys
{

    class VDWForce
    : public Force
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        boost::shared_ptr< mol::SimpleMolecule< mol::Atom> >      m_Molecule;
        float                                                     m_Cutoff;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        VDWForce()
        : Force(),
        m_Y()
        {}

        //! construct from data
        VDWForce
        (
                const float WEIGHT,
                const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE,
                const float &CUTOFF = std::numeric_limits< float>::max()
        )
        : Force( WEIGHT),
        m_Molecule( IMPLICIT_MOLECULE),
        m_Cutoff( CUTOFF)
        {}

        //! copy constructor
        VDWForce( const ElectrostaticForce &ORIGINAL)
        : Force( ORIGINAL),
        m_Molecule( ORIGINAL.m_Molecule),
        m_Cutoff( ORIGINAL.m_Cutoff)
        {}

        //! virtual destructor
        virtual ~VDWForce(){}

        //! virtual copy constructor
        virtual VDWForce *Clone() const{ return new VDWForce( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////
        math::Vector3N operator()(  const mol::Atom &ATOM_A, const mol::Atom &ATOM_B) const
        {
            math::Vector3N connector( ATOM_B->GetPosition() - ATOM_A->GetPosition());
            float distance( connector.Length());
            if( distance > m_Cutoff || ATOM_A->GetType() == "H" || ATOM_B->GetType() == "H")
            { return math::Vector3N( new math::Vector3N( 0.0, 0.0, 0.0));}
            return
                math::Vector3N
                (
                        new math::Vector3N( Force::m_Weight * 12.0 * sqrt( ATOM_A->GetEpsilon()) * ATOM_B->GetPartialCharge() * std::pow( distance, -3) * connector)
                );
        }


        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////
        boost::shared_ptr< phys::PotentialForce> ElectrostaticForce::PotentialForce( const math::Vector3N & POS) const
        {
    //        std::cout << __PRETTY_FUNCTION__ << std::endl;
            math::Vector3N potential( new math::Vector3N());
            for( std::vector< mol::Atom >::const_iterator mol_itr = m_Molecule->GetAtoms().begin(); mol_itr != m_Molecule->GetAtoms().end(); ++mol_itr)
            {
                math::Vector3N connector( ( *mol_itr)->GetPosition() - *POS);
                float distance( connector.GetLengthAndNormalize());
                std::cout << "type: " << ( *mol_itr)->GetType() << " partial_charge: " << ( *mol_itr)->GetPartialCharge() << " distance: " << distance << " <std::pow( distance, -3) * 332.0716 * partial_charge: " << 332.0716 * ( *mol_itr)->GetPartialCharge() * std::pow( distance, -3) << "> (" << std::pow( distance, -3) << ")" << std::endl;
                if( distance <= m_Cutoff || ( *mol_itr)->GetType() != "H")  // TODO !! why no H??
                {
                    *potential += 332.0716 * ( *mol_itr)->GetPartialCharge() * std::pow( distance, -3) * connector;  // TODO !!
    //                *potential += Force::m_Weight * 332.0716  * std::pow( distance, -3) * connector.NormalizedCopy();
                }
            }
            std::cout << std::endl << std::endl;
    //        std::cout << "potential: " << potential->Length() << " times " << Force::m_Weight << " = ";
            *potential *= Force::m_Weight;
    //        std::cout << "potential: " << potential->Length() << std::endl;
            return boost::shared_ptr< phys::PotentialForce>( new PotentialElectrostaticForce( potential));
        }

        boost::shared_ptr< phys::PotentialForce> ElectrostaticForce::GetAssociatedPotentialForceObject() const
        { /*std::cout << __PRETTY_FUNCTION__ << std::endl;*/ return boost::shared_ptr< phys::PotentialForce>( new PotentialElectrostaticForce());}

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

    }; // end class VDWForce
} // end namespace phys




#endif /* VDW_ENERGY_H_ */

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


#ifndef ELECTROSTATIC_FORCE_H_
#define ELECTROSTATIC_FORCE_H_


#include "force.h"
#include "potential_electrostatic_force.h"
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"


namespace phys
{

    class ElectrostaticForce
    : public Force
    {
    protected:
    	static float										  sm_ForceEnergyUnitFactor;
        boost::shared_ptr< mol::SimpleMolecule< mol::Atom> >  m_Molecule;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        ElectrostaticForce( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE, const float &CUTOFF = std::numeric_limits< float>::max());

        ElectrostaticForce( const float &CUTOFF = std::numeric_limits< float>::max());

        //! copy constructor
        ElectrostaticForce( const ElectrostaticForce &ORIGINAL);

        //! virtual destructor
        virtual ~ElectrostaticForce();

        //! virtual copy constructor
        virtual ElectrostaticForce *Clone() const;

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        //! the force caused by ATOM_B acting on ATOM_A
        virtual math::Vector3N operator()( const mol::Atom &ATOM_A, const mol::Atom &ATOM_B) const;

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual boost::shared_ptr< phys::PotentialForce> PotentialForce( const math::Vector3N & POS) const;

        virtual boost::shared_ptr< phys::PotentialForce> GetAssociatedPotentialForceObject() const;

        virtual void SetKjouleNmUnits();

        /////////////////////////
        //      Read/Write     //
        /////////////////////////


        virtual std::istream &Read( std::istream &STREAM);

        virtual std::ostream &Write( std::ostream &STREAM) const;

        virtual std::string GetClassName() const;

    }; // end class ElectrostaticForce
} // end namespace phys




#endif /* ELECTROSTATIC_FORCE_H_ */

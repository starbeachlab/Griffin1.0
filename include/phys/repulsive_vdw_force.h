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


#ifndef REPULSIVE_VDW_FORCE_H_
#define REPULSIVE_VDW_FORCE_H_


#include "force.h"
#include "potential_repulsive_vdw_force.h"
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"


namespace phys
{

    class RepulsiveVanDerWaalsForce
    : public Force
    {

    protected:
//    	static float										  sm_ForceUnitFactor;
//    	static float    									  sm_EnergyUnitFactor;
        boost::shared_ptr< mol::SimpleMolecule< mol::Atom> >  m_Molecule;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        RepulsiveVanDerWaalsForce( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE, const float &CUTOFF = std::numeric_limits< float>::max());

        //! default constructor
        RepulsiveVanDerWaalsForce( const float &CUTOFF = std::numeric_limits< float>::max());

        //! copy constructor (meaningless)
        RepulsiveVanDerWaalsForce( const RepulsiveVanDerWaalsForce &ORIGINAL);

        //! virtual destructor
        virtual ~RepulsiveVanDerWaalsForce();

        //! virtual copy constructor
        virtual RepulsiveVanDerWaalsForce *Clone() const;

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////
        virtual math::Vector3N operator()( const mol::Atom &ATOM_A, const mol::Atom &ATOM_B) const;

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual boost::shared_ptr< phys::PotentialForce>  PotentialForce( const math::Vector3N & POS) const;

        //virtual math::Vector3N PotentialForce( const math::Vector3N & POS, const boost::shared_ptr< mol::MoleculeIF> &MOLECULE);

        virtual boost::shared_ptr< phys::PotentialForce> GetAssociatedPotentialForceObject() const;

//        virtual void SetKjouleNmUnits();

        /////////////////////////
        //      Read/Write     //
        /////////////////////////
        virtual std::istream &Read( std::istream &STREAM);

        virtual std::ostream &Write( std::ostream &STREAM) const;

        virtual std::string GetClassName() const;

    }; // end class RepulsiveVanDerWaalsForce
} // end namespace phys




#endif /* REPULSIVE_VDW_FORCE_H_*/

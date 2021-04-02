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


#ifndef POTENTIAL_attractive_vdw_force_H_
#define POTENTIAL_attractive_vdw_force_H_


#include "../readwrite/class_name.h"
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"
#include "../math/vector3N.h"
#include "../molecule/atom.h"
#include "potential_force.h"

namespace phys
{

    class PotentialAttractiveVanDerWaalsForce
    : public PotentialForce
    {
    private:
    	static math::Vector3N   s_Local;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        PotentialAttractiveVanDerWaalsForce();

        //! construct from data
        PotentialAttractiveVanDerWaalsForce( const math::Vector3N &FORCE, const float &ENERGY);

        //! copy constructor
        PotentialAttractiveVanDerWaalsForce( const PotentialAttractiveVanDerWaalsForce &ORIGINAL);

        //! virtual destructor
        virtual ~PotentialAttractiveVanDerWaalsForce();

        //! virtual copy constructor
        virtual PotentialAttractiveVanDerWaalsForce *Clone() const;

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        virtual math::Vector3N operator() ( const mol::Atom &ATOM) const;


        virtual float Energy( const mol::Atom &ATOM) const;


        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM);

        virtual std::ostream &Write( std::ostream &STREAM) const;

        virtual std::string GetClassName() const;

    }; // end class PotentialAttractiveVanDerWaalsForce


} // end namespace phys




#endif /* VanDerWaals_FORCE_H_ */

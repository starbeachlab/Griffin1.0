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


#ifndef FORCE_FACTORY_H_
#define FORCE_FACTORY_H_

#include "potential_force_container.h"
#include "force_container.h"
#include "electrostatic_force.h"
#include "attractive_vdw_force.h"
#include "repulsive_vdw_force.h"
//#include "surf_force.h"
//#include "density_force.h"
#include "../string/string_functions.h"
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"


namespace phys
{

    namespace factory
    {
        std::istream &ReadForces( std::istream &STREAM, boost::shared_ptr< phys::ForceContainer> &FORCES, const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE);

        boost::shared_ptr< PotentialForceContainer>
        BuildPotentialForceContainerFromForceContainer( const phys::ForceContainer &FORCES);


        PotentialForceContainer
        BuildPotentialForces
        (
                const boost::shared_ptr< phys::ForceContainer> &FORCES,
                const math::Vector3N &GRID_POS,
                const float &MIN_FORCE_MAGNITUDE,
                const float &DISTANCE
//                const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOL
        );

    } // end namespace factory
} // end namespace phys




#endif /* FORCE_FACTORY_H_ */

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


#ifndef MOLECULE_FORCE_GRID_FACTORY_H_
#define MOLECULE_FORCE_GRID_FACTORY_H_

#include "molecule_force_grid.h"
#include "../storage/interpol_grid.h"

namespace mol
{
    namespace factory
    {
		void
		RescaleConstantForces
		(
				boost::shared_ptr< mol::MoleculeForceGrid> &GRID,
				const float &VALUE
		);

		void
		RescaleForces
		(
				boost::shared_ptr< mol::MoleculeForceGrid> &GRID,
				const std::string &FORCE_TYPE,
//				const std::string &SCALING_STYLE,
				const float &VALUE
		);

		void
		BuildGridFromParallelFiles
		(
				boost::shared_ptr< mol::MoleculeForceGrid> &GRID,
				const CommandLineManager &CMD

		);

		boost::shared_ptr< mol::MoleculeForceGrid>
		MergeGrids
		(
				const std::vector< boost::shared_ptr< mol::MoleculeForceGrid> > &SUB_GRIDS,
				const CommandLineManager &CMD
		);

		void
		ReplaceInteractionWithTransientGridPoints
		(
				boost::shared_ptr< mol::MoleculeForceGrid> &GRID
		);


    } // end namespace factory
} // end namespace mol


#endif /* MOLECULE_FORCE_GRID_FACTORY_H_ */

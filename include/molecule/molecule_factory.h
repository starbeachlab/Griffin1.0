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


#ifndef MOLECULE_FACTORY_H_
#define MOLECULE_FACTORY_H_

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"
#include "../storage/shared_pointer_vector.t.h"
#include "../storage/vector3N.t.h"
#include "simple_molecule.t.h"
#include "atom.h"
#include "simple_molecule_file_handler.t.h"
#include "../command/command_line_manager.h"
#include "../readwrite/stream_functions.h"
#include "../storage/map.t.h"


namespace mol
{
    namespace factory
    {
		void Read( const store::ShPtrVec< SimpleMolecule< Atom> > &MOL, const std::istream &READ, const std::string &INPUT_FILE_FORMAT);

		void CalculateGridLimits( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOL, math::Vector3N &MIN, math::Vector3N &MAX, const math::Vector3N &DELTA, const float &CUTOFF);

		store::Vector3N< size_t> CalculateNumberOfBins( const math::Vector3N &MIN, const math::Vector3N &MAX, const math::Vector3N &DELTA);

		void AdjustLimits( const float &PRECISION, const math::Vector3N &MIN, math::Vector3N &MAX, const math::Vector3N &DELTA);

		void AdjustLimits( math::Vector3N &MIN, math::Vector3N &MAX, math::Vector3N &DELTA, const math::Vector3N &FIXED_MIN, const math::Vector3N &FIXED_MAX, const math::Vector3N &FIXED_DELTA);

		void
		BuildImplicitMolecule
		(
				const CommandLineManager &COMMAND,
				boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE
		);


		void
		CalcVdwFactors( boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE);

    } // end namespace factory
} // end namespace mol




#endif /* MOLECULE_FACTORY_H_ */

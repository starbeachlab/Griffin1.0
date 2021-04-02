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


#ifndef ENUM_HANDLER_ENERGY_H
#define ENUM_HANDLER_ENERGY_H

#include "energy_enum.h"
#include "enum_handler.t.h"

namespace util
{
	template<>
	void
	BuildEnumHandler< EnergyEnum>()
	{
	    StandardWrite( __PRETTY_FUNCTION__);
		EnumHandler< EnergyEnum> hand;
		hand.Insert( e_AttractiveVdw, "vdw-attractive");
		hand.Insert( e_RepulsiveVdw, "vdw-repulsive");
		hand.Insert( e_Coulomb, "coulomb");
		hand.Insert( e_Density, "density");
		hand.Insert( e_UNDEFINED_Energy,  "util::UNDEFINED_Energy");
	}
} // end namespace util

#endif


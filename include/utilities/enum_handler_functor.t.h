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


#ifndef ENUM_HANDLER_FUNCTOR_H
#define ENUM_HANDLER_FUNCTOR_H

#include "enum_handler.t.h"
#include "functor_enum.h"

namespace util
{
	template<>
	void
	BuildEnumHandler< FunctorEnum>()
	{
	    StandardWrite( __PRETTY_FUNCTION__);
		EnumHandler< FunctorEnum> hand;
		hand.Insert( e_Function, "math::Function");
		hand.Insert( e_VoidGridPoint, "store::VoidGridPoint");
		hand.Insert( e_SurfGridPoint, "store::SurfGridPoint");
		hand.Insert( e_ConstantGridPoint, "store::ConstantGridPoint");
		hand.Insert( e_InteractionGridPoint, "store::InteractionGridPoint");
		hand.Insert( e_TypeMappedGridPoint, "store::TypeMappedGridPoint");
		hand.Insert( e_RecursiveTypeMappedGridPoint, "store::RecursiveTypeMappedGridPoint");
		hand.Insert( e_SumGridPoint, "store::SumGridPoint");
		hand.Insert( e_TransientInteractionGridPoint, "store::TransientInteractionGridPoint");
		hand.Insert( e_PotentialElectrostaticForce, "phys::PotentialElectrostaticForce");
		hand.Insert( e_PotentialAttractiveVanDerWaalsForce, "phys::PotentialAttractiveVanDerWaalsForce");
		hand.Insert( e_PotentialRepulsiveVanDerWaalsForce, "phys::PotentialRepulsiveVanDerWaalsForce");
		hand.Insert( e_UNDEFINED_Functor,  "util::UNDEFINED_Functor");
	}
} // end namespace util

#endif

